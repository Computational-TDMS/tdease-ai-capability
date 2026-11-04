//! Output envelope and error JSON shapes for the tdeasenexus-toolbox contract.
//!
//! Success: one object on stdout (`{schemaVersion, tool, version, command,
//! params, source, elapsedMs, warnings, result}`). Errors: one object on
//! stderr (`{error: {code, message, details?}}`), exit 1 for runtime errors,
//! exit 2 for usage errors. stdout stays empty on error — result and error
//! never mix.

use serde::Serialize;
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

/// The standard success envelope. `source` is absent when the command read no
/// files; `elapsed_ms` is filled in by the caller after dispatch.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct OutputEnvelope {
    pub schema_version: u32,
    pub tool: &'static str,
    pub version: &'static str,
    pub command: String,
    pub params: Value,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub source: Vec<SourceInfo>,
    pub elapsed_ms: u128,
    pub warnings: Vec<String>,
    pub result: Value,
}

/// One input file the command actually read.
#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SourceInfo {
    pub path: String,
    pub size: u64,
    pub digest: SourceDigest,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct SourceDigest {
    pub algo: String,
    pub value: String,
}

impl OutputEnvelope {
    pub fn new(command: &str, params: Value, source: Vec<SourceInfo>, result: Value) -> Self {
        Self {
            schema_version: 1,
            tool: "tdeasenexus-toolbox",
            version: env!("CARGO_PKG_VERSION"),
            command: command.to_owned(),
            params,
            source,
            elapsed_ms: 0,
            warnings: Vec::new(),
            result,
        }
    }
}

/// sha256 over the first 1 MiB (fast provenance marker — labelled as such).
pub fn head_digest(bytes: &[u8]) -> String {
    let head = &bytes[..bytes.len().min(1 << 20)];
    hex(&Sha256::digest(head))
}

/// sha256 over the whole file.
pub fn full_digest(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

fn hex(digest: &[u8]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

pub fn emit_ok(envelope: &OutputEnvelope) {
    println!("{}", serde_json::to_string(envelope).unwrap());
}

/// Emit a runtime error as stderr JSON, exit 1.
pub fn emit_json_error(error: &Value) {
    eprintln!("{}", serde_json::to_string(error).unwrap());
}

/// Emit a usage error (bad arguments) as stderr JSON, exit 2.
pub fn emit_usage_error(message: &str) {
    emit_json_error(&json!({ "error": { "code": "USAGE", "message": message } }));
}
