//! A minimal **mzTab-M** writer (PSI standard, v2.0.0-M).
//!
//! mzTab is the HUPO-PSI tab-separated exchange format for processed
//! proteomics results. This module emits the subset that carries our
//! fragmentation / candidate-discrimination results as an external-readable
//! intermediate file: an `MTD` metadata header plus the **small-molecule**
//! table (`SMH` header + `SML` rows — the mzTab 2.0 molecular-table naming
//! that `pyteomics.mztab` parses; the mzTab 1.0 names for the same table were
//! `MZH`/`MZM`).
//!
//! The writer is deliberately generic — it formats key/value metadata and a
//! rectangular table; the caller decides which columns carry the result
//! (sequence, modifications, theoretical mass, matched fragments, …).
//! It has no dependency on the rest of the crate beyond the std library.

/// Escape a single cell: mzTab values are tab-separated, so a literal tab,
/// CR or LF inside a value is replaced with a space (the PSI spec has no
/// escaping mechanism; it just forbids the separator in values).
fn escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '\t' | '\r' | '\n' => out.push(' '),
            other => out.push(other),
        }
    }
    out
}

/// Builds a well-formed mzTab document line by line.
#[derive(Debug, Clone, Default)]
pub struct MzTabWriter {
    metas: Vec<(String, String)>,
    molecular_columns: Vec<String>,
    molecular_rows: Vec<Vec<String>>,
}

impl MzTabWriter {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add an `MTD <key>\t<value>` metadata line.
    pub fn meta(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.metas.push((key.into(), value.into()));
    }

    /// Define the molecular-table columns (the `MZH` header row). Every
    /// `add_molecular_row` must supply the same number of cells.
    pub fn set_molecular_columns(&mut self, columns: Vec<impl Into<String>>) {
        self.molecular_columns = columns.into_iter().map(Into::into).collect();
    }

    /// Append one molecular row (`MZM`). Cell count must match the columns.
    pub fn add_molecular_row(&mut self, row: Vec<impl Into<String>>) {
        debug_assert_eq!(
            row.len(),
            self.molecular_columns.len(),
            "mzTab MZM row width must match MZH header width"
        );
        self.molecular_rows
            .push(row.into_iter().map(Into::into).collect());
    }

    /// Render the complete tab-separated mzTab document.
    pub fn render(&self) -> String {
        let mut out = String::new();
        for (key, value) in &self.metas {
            out.push_str("MTD\t");
            out.push_str(&escape(key));
            out.push('\t');
            out.push_str(&escape(value));
            out.push('\n');
        }
        if self.molecular_columns.is_empty() {
            return out;
        }
        // mzTab 2.0 small-molecule table naming (pyteomics-compatible).
        out.push_str("SMH\t");
        out.push_str(
            &self
                .molecular_columns
                .iter()
                .map(|c| escape(c))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        out.push('\n');
        for row in &self.molecular_rows {
            out.push_str("SML\t");
            out.push_str(&row.iter().map(|c| escape(c)).collect::<Vec<_>>().join("\t"));
            out.push('\n');
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_standard_mztab_document() {
        let mut mz = MzTabWriter::new();
        mz.meta("mzTab-version", "2.0.0-M");
        mz.meta("mzTab-mode", "Complete");
        mz.meta("mzTab-type", "Identification");
        mz.meta("description", "fragment-theoretical");
        mz.set_molecular_columns(vec!["mzTab_id", "chemical_formula", "opt_mzTab_sequence"]);
        mz.add_molecular_row(vec!["1", "C9H8NO2Br", "SGRG"]);
        let doc = mz.render();
        assert!(doc.starts_with("MTD\tmzTab-version\t2.0.0-M\n"));
        assert!(doc.contains("SMH\tmzTab_id\tchemical_formula\topt_mzTab_sequence\n"));
        assert!(doc.contains("SML\t1\tC9H8NO2Br\tSGRG\n"));
    }

    #[test]
    fn escapes_separator_chars_in_values() {
        let mut mz = MzTabWriter::new();
        mz.set_molecular_columns(vec!["a", "b"]);
        mz.add_molecular_row(vec!["x\ty", "z\nw"]);
        let doc = mz.render();
        assert!(doc.contains("SML\tx y\tz w\n"), "got: {doc}");
    }
}
