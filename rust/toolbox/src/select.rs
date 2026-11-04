use mzdata::io::SpectrumSource;
use mzdata::params::{ParamDescribed, ParamValue};
use mzdata::prelude::SpectrumLike;
use serde::{Deserialize, Serialize};

use crate::mzml::{open_reader, scan_token};
use crate::{MassSpecError, Peak};

/// The mzML header's total ion current (MS:1000285), when the file declares it.
pub(crate) fn header_total_ion_current(params: &[mzdata::params::Param]) -> Option<f64> {
    params
        .iter()
        .find(|param| param.name == "total ion current")
        .and_then(|param| param.to_f64().ok())
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase", tag = "field")]
pub enum SelectPredicate {
    MsLevel {
        value: u8,
    },
    RetentionTime {
        min: f64,
        max: f64,
    },
    PrecursorMz {
        min: f64,
        max: f64,
    },
    /// Precursor neutral mass (selected m/z × z − z·proton, from the scan's
    /// own charge header). Unlike `PrecursorMz` this is one interval across
    /// charges — the top-down family window the original IsoForma
    /// `pull_scan_numbers` filtered on at the metadata stage.
    PrecursorNeutralMass {
        min: f64,
        max: f64,
    },
    Charge {
        value: i32,
    },
    Activation {
        value: String,
    },
    ProductMz {
        min: f64,
        max: f64,
    },
    NeutralLoss {
        min: f64,
        max: f64,
    },
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SelectPlan {
    pub predicates: Vec<SelectPredicate>,
    pub limit: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrecursorMetadata {
    pub selected_mz: f64,
    pub intensity: f32,
    pub charge: Option<i32>,
    pub isolation_target: f32,
    pub isolation_lower_bound: f32,
    pub isolation_upper_bound: f32,
    pub activation: Option<String>,
    #[serde(default)]
    pub fragment_ion_types: Option<Vec<crate::fragment::IonType>>,
    pub collision_energy: Option<f32>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SelectedSpectrum {
    pub id: String,
    pub scan: String,
    pub index: usize,
    pub ms_level: u8,
    pub retention_time: f64,
    pub retention_time_unit: String,
    /// The spectrum's total ion current: the mzML header value
    /// (`total ion current`, MS:1000285), falling back to the summed peak
    /// intensities when the header omits it.
    pub total_ion_current: f64,
    /// `None` when peaks were not decoded — statements without a
    /// `product_mz`/`neutral_loss` predicate don't need peak-level evidence.
    pub peak_count: Option<usize>,
    pub centroided: bool,
    pub precursor: Option<PrecursorMetadata>,
    pub matched_products: Vec<Peak>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct SelectQueryResult {
    pub plan: SelectPlan,
    pub scanned_spectrum_count: usize,
    pub matched_spectrum_count: usize,
    pub truncated: bool,
    pub spectra: Vec<SelectedSpectrum>,
}

pub fn query_select_mzml(
    bytes: &[u8],
    statement: &str,
) -> Result<SelectQueryResult, MassSpecError> {
    let plan = parse_select(statement)?;
    let needs_peak_evidence = plan.predicates.iter().any(|predicate| {
        matches!(
            predicate,
            SelectPredicate::ProductMz { .. } | SelectPredicate::NeutralLoss { .. }
        )
    });
    let mut reader = open_reader(bytes)?;
    let mut scanned_spectrum_count = 0;
    let mut spectra = Vec::new();
    let mut truncated = false;

    for mut spectrum in reader.iter() {
        scanned_spectrum_count += 1;
        let precursor = spectrum.precursor().and_then(|value| {
            let ion = value.ions.first()?;
            let activation = value
                .activation
                .method()
                .map(|method| method.name().to_owned());
            let fragment_ion_types =
                crate::fragment::ion_types_for_activation(activation.as_deref()).ok();
            Some(PrecursorMetadata {
                selected_mz: ion.mz,
                intensity: ion.intensity,
                charge: ion.charge,
                isolation_target: value.isolation_window.target,
                isolation_lower_bound: value.isolation_window.lower_bound,
                isolation_upper_bound: value.isolation_window.upper_bound,
                activation,
                fragment_ion_types,
                collision_energy: (value.activation.energy != 0.0)
                    .then_some(value.activation.energy),
            })
        });
        let retention_time = spectrum.start_time();
        if !metadata_matches(
            &plan,
            spectrum.ms_level(),
            retention_time,
            precursor.as_ref(),
        ) {
            continue;
        }

        // Peak decoding dominates runtime; run it only when the statement
        // needs peak-level evidence or the header omits the total ion current.
        let header_tic = header_total_ion_current(spectrum.params());
        let peaks_built = needs_peak_evidence || header_tic.is_none();
        if peaks_built {
            spectrum
                .try_build_peaks()
                .map_err(|error| MassSpecError::PeakDecode(error.to_string()))?;
        }
        let total_ion_current = match header_tic {
            Some(value) => value,
            None => spectrum
                .peaks()
                .iter()
                .map(|peak| peak.intensity as f64)
                .sum(),
        };
        let (peak_count, matched_products) = if peaks_built {
            let peak_data = spectrum.peaks();
            let centroid_peaks: Vec<Peak> = peak_data
                .iter()
                .map(|peak| Peak {
                    mz: peak.mz,
                    intensity: peak.intensity,
                })
                .collect();
            let matched = match peak_matches(&plan, precursor.as_ref(), &centroid_peaks) {
                Some(peaks) => peaks,
                None => continue,
            };
            (Some(peak_data.len()), matched)
        } else {
            (None, Vec::new())
        };
        if spectra.len() == plan.limit {
            truncated = true;
            break;
        }
        spectra.push(SelectedSpectrum {
            id: spectrum.id().to_owned(),
            scan: scan_token(spectrum.id()).unwrap_or_else(|| spectrum.id().to_owned()),
            index: spectrum.index(),
            ms_level: spectrum.ms_level(),
            retention_time,
            retention_time_unit: "minute".to_owned(),
            total_ion_current,
            peak_count,
            centroided: spectrum.signal_continuity()
                == mzdata::spectrum::SignalContinuity::Centroid,
            precursor,
            matched_products,
        });
    }

    Ok(SelectQueryResult {
        plan,
        scanned_spectrum_count,
        matched_spectrum_count: spectra.len(),
        truncated,
        spectra,
    })
}

fn metadata_matches(
    plan: &SelectPlan,
    ms_level: u8,
    retention_time: f64,
    precursor: Option<&PrecursorMetadata>,
) -> bool {
    plan.predicates.iter().all(|predicate| match predicate {
        SelectPredicate::MsLevel { value } => ms_level == *value,
        SelectPredicate::RetentionTime { min, max } => {
            retention_time.is_finite() && retention_time >= *min && retention_time <= *max
        }
        SelectPredicate::PrecursorMz { min, max } => {
            precursor.is_some_and(|value| value.selected_mz >= *min && value.selected_mz <= *max)
        }
        SelectPredicate::PrecursorNeutralMass { min, max } => {
            precursor.is_some_and(|value| match value.charge {
                Some(charge) => {
                    let z = charge.abs() as f64;
                    let neutral = value.selected_mz * z - z * crate::constants::PROTON;
                    neutral >= *min && neutral <= *max
                }
                None => false,
            })
        }
        SelectPredicate::Charge { value } => {
            precursor.is_some_and(|precursor| precursor.charge == Some(*value))
        }
        SelectPredicate::Activation { value } => precursor
            .and_then(|precursor| precursor.activation.as_deref())
            .is_some_and(|activation| activation.eq_ignore_ascii_case(value)),
        SelectPredicate::ProductMz { .. } | SelectPredicate::NeutralLoss { .. } => true,
    })
}

fn peak_matches(
    plan: &SelectPlan,
    precursor: Option<&PrecursorMetadata>,
    peaks: &[Peak],
) -> Option<Vec<Peak>> {
    let mut matched = Vec::new();
    for predicate in &plan.predicates {
        let range = match predicate {
            SelectPredicate::ProductMz { min, max } => (*min, *max),
            SelectPredicate::NeutralLoss { min, max } => {
                let precursor_mz = precursor?.selected_mz;
                (precursor_mz - *max, precursor_mz - *min)
            }
            _ => continue,
        };
        let evidence: Vec<Peak> = peaks
            .iter()
            .filter(|peak| peak.mz >= range.0 && peak.mz <= range.1)
            .map(|peak| Peak {
                mz: peak.mz,
                intensity: peak.intensity,
            })
            .collect();
        if evidence.is_empty() {
            return None;
        }
        for peak in evidence {
            if !matched.iter().any(|existing: &Peak| existing.mz == peak.mz) {
                matched.push(peak);
            }
        }
    }
    matched.sort_unstable_by(|left, right| left.mz.total_cmp(&right.mz));
    Some(matched)
}

#[derive(Debug, Clone, PartialEq)]
enum Token {
    Word(String),
    Number(f64),
    Text(String),
    Eq,
    Semicolon,
}

fn parse_select(statement: &str) -> Result<SelectPlan, MassSpecError> {
    let tokens = lex(statement)?;
    let mut parser = Parser {
        tokens,
        position: 0,
    };
    parser.expect_word("select")?;
    parser.expect_word("spectra")?;
    parser.expect_word("from")?;
    parser.expect_word("mzml")?;
    parser.expect_word("where")?;
    let mut predicates = Vec::new();
    loop {
        predicates.push(parser.predicate()?);
        if parser.consume_word("and") {
            continue;
        }
        break;
    }
    let limit = if parser.consume_word("limit") {
        let value = parser.number()?;
        if value.fract() != 0.0 || !(1.0..=10_000.0).contains(&value) {
            return invalid("LIMIT must be an integer between 1 and 10000");
        }
        value as usize
    } else {
        100
    };
    parser.consume_semicolon();
    if parser.position != parser.tokens.len() {
        return invalid("unexpected tokens after SELECT statement");
    }
    Ok(SelectPlan { predicates, limit })
}

struct Parser {
    tokens: Vec<Token>,
    position: usize,
}

impl Parser {
    fn predicate(&mut self) -> Result<SelectPredicate, MassSpecError> {
        let field = self.word()?;
        match field.as_str() {
            "ms_level" => {
                self.expect_eq()?;
                let value = self.number()?;
                if value.fract() != 0.0 || !(1.0..=255.0).contains(&value) {
                    return invalid("ms_level must be an integer between 1 and 255");
                }
                Ok(SelectPredicate::MsLevel { value: value as u8 })
            }
            "charge" => {
                self.expect_eq()?;
                let value = self.number()?;
                if value.fract() != 0.0 || value < i32::MIN as f64 || value > i32::MAX as f64 {
                    return invalid("charge must be an integer");
                }
                Ok(SelectPredicate::Charge {
                    value: value as i32,
                })
            }
            "activation" => {
                self.expect_eq()?;
                Ok(SelectPredicate::Activation {
                    value: self.text()?,
                })
            }
            "rt" | "retention_time" => {
                let (min, max) = self.range()?;
                Ok(SelectPredicate::RetentionTime { min, max })
            }
            "precursor_mz" => {
                let (min, max) = self.range()?;
                Ok(SelectPredicate::PrecursorMz { min, max })
            }
            "precursor_neutral_mass" => {
                let (min, max) = self.range()?;
                Ok(SelectPredicate::PrecursorNeutralMass { min, max })
            }
            "product_mz" => {
                let (min, max) = self.range()?;
                Ok(SelectPredicate::ProductMz { min, max })
            }
            "neutral_loss" => {
                let (min, max) = self.range()?;
                Ok(SelectPredicate::NeutralLoss { min, max })
            }
            _ => invalid(&format!("unsupported spectra field: {field}")),
        }
    }

    fn range(&mut self) -> Result<(f64, f64), MassSpecError> {
        if self.consume_word("within") {
            let center = self.number()?;
            let unit = self.word()?;
            let tolerance = self.number()?;
            if !center.is_finite() || center <= 0.0 || !tolerance.is_finite() || tolerance <= 0.0 {
                return invalid("WITHIN requires a positive center and tolerance");
            }
            let width = match unit.as_str() {
                "da" => tolerance,
                "ppm" => center * tolerance / 1_000_000.0,
                _ => return invalid("WITHIN unit must be DA or PPM"),
            };
            return Ok((center - width, center + width));
        }
        self.expect_word("between")?;
        let min = self.number()?;
        self.expect_word("and")?;
        let max = self.number()?;
        if !min.is_finite() || !max.is_finite() || min > max {
            return invalid("BETWEEN requires an ordered finite range");
        }
        Ok((min, max))
    }

    fn word(&mut self) -> Result<String, MassSpecError> {
        match self.tokens.get(self.position) {
            Some(Token::Word(value)) => {
                self.position += 1;
                Ok(value.clone())
            }
            _ => invalid("expected an identifier"),
        }
    }

    fn text(&mut self) -> Result<String, MassSpecError> {
        match self.tokens.get(self.position) {
            Some(Token::Text(value)) => {
                self.position += 1;
                Ok(value.clone())
            }
            _ => invalid("expected a quoted string"),
        }
    }

    fn number(&mut self) -> Result<f64, MassSpecError> {
        match self.tokens.get(self.position) {
            Some(Token::Number(value)) => {
                self.position += 1;
                Ok(*value)
            }
            _ => invalid("expected a number"),
        }
    }

    fn expect_word(&mut self, expected: &str) -> Result<(), MassSpecError> {
        if self.consume_word(expected) {
            Ok(())
        } else {
            invalid(&format!("expected {expected}"))
        }
    }

    fn consume_word(&mut self, expected: &str) -> bool {
        match self.tokens.get(self.position) {
            Some(Token::Word(value)) if value == expected => {
                self.position += 1;
                true
            }
            _ => false,
        }
    }

    fn expect_eq(&mut self) -> Result<(), MassSpecError> {
        if matches!(self.tokens.get(self.position), Some(Token::Eq)) {
            self.position += 1;
            Ok(())
        } else {
            invalid("expected =")
        }
    }

    fn consume_semicolon(&mut self) {
        if matches!(self.tokens.get(self.position), Some(Token::Semicolon)) {
            self.position += 1;
        }
    }
}

fn lex(input: &str) -> Result<Vec<Token>, MassSpecError> {
    let chars: Vec<char> = input.chars().collect();
    let mut tokens = Vec::new();
    let mut index = 0;
    while index < chars.len() {
        match chars[index] {
            value if value.is_whitespace() => index += 1,
            '=' => {
                tokens.push(Token::Eq);
                index += 1;
            }
            ';' => {
                tokens.push(Token::Semicolon);
                index += 1;
            }
            '\'' => {
                index += 1;
                let mut value = String::new();
                loop {
                    match chars.get(index) {
                        Some('\'') if chars.get(index + 1) == Some(&'\'') => {
                            value.push('\'');
                            index += 2;
                        }
                        Some('\'') => {
                            index += 1;
                            break;
                        }
                        Some(character) => {
                            value.push(*character);
                            index += 1;
                        }
                        None => return invalid("unterminated quoted string"),
                    }
                }
                tokens.push(Token::Text(value));
            }
            value if value.is_ascii_digit() || matches!(value, '+' | '-' | '.') => {
                let start = index;
                index += 1;
                while index < chars.len()
                    && (chars[index].is_ascii_digit()
                        || matches!(chars[index], '.' | 'e' | 'E' | '+' | '-'))
                {
                    index += 1;
                }
                let text: String = chars[start..index].iter().collect();
                let value = text
                    .parse::<f64>()
                    .map_err(|_| MassSpecError::InvalidQuery(format!("invalid number: {text}")))?;
                tokens.push(Token::Number(value));
            }
            value if value.is_ascii_alphabetic() || value == '_' => {
                let start = index;
                index += 1;
                while index < chars.len()
                    && (chars[index].is_ascii_alphanumeric() || chars[index] == '_')
                {
                    index += 1;
                }
                tokens.push(Token::Word(
                    chars[start..index]
                        .iter()
                        .collect::<String>()
                        .to_ascii_lowercase(),
                ));
            }
            value => return invalid(&format!("unexpected character: {value}")),
        }
    }
    Ok(tokens)
}

fn invalid<T>(message: &str) -> Result<T, MassSpecError> {
    Err(MassSpecError::InvalidQuery(message.to_owned()))
}
