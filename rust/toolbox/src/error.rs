use thiserror::Error;

#[derive(Debug, Error)]
pub enum MassSpecError {
    #[error("failed to open mass-spec data: {0}")]
    Open(String),
    #[error("no spectrum matched {0}")]
    SpectrumNotFound(String),
    #[error("failed to decode spectrum peaks: {0}")]
    PeakDecode(String),
    #[error("invalid scan query: {0}")]
    InvalidQuery(String),
}

impl MassSpecError {
    pub fn code(&self) -> &'static str {
        match self {
            Self::Open(_) => "MASS_SPEC_OPEN_FAILED",
            Self::SpectrumNotFound(_) => "MASS_SPEC_SPECTRUM_NOT_FOUND",
            Self::PeakDecode(_) => "MASS_SPEC_PEAK_DECODE_FAILED",
            Self::InvalidQuery(_) => "MASS_SPEC_INVALID_QUERY",
        }
    }
}
