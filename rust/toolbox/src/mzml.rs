use std::io::Cursor;

use mzdata::io::MZReaderType;
use mzpeaks::{CentroidPeak, DeconvolutedPeak};

use crate::MassSpecError;

pub type Reader = MZReaderType<Cursor<Vec<u8>>, CentroidPeak, DeconvolutedPeak>;

pub fn open_reader(bytes: &[u8]) -> Result<Reader, MassSpecError> {
    Reader::open_read_seek(Cursor::new(bytes.to_vec()))
        .map_err(|error| MassSpecError::Open(error.to_string()))
}

pub(crate) fn scan_token(id: &str) -> Option<String> {
    let tail = id.rsplit_once("scan=")?.1;
    let token: String = tail.chars().take_while(char::is_ascii_digit).collect();
    (!token.is_empty()).then_some(token)
}
