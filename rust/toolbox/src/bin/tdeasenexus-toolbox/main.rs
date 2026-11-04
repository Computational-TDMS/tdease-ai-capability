//! tdeasenexus-toolbox — the native adapter for the mass-spec core.
//!
//! The single delivery mechanism for the `rust/toolbox` core: one JSON
//! object per invocation on stdout; stable `code`/`message` errors on stderr;
//! `--help` is the interface. Internal agents run this via bash and Python
//! skills call it via `subprocess`; neither re-implements mzML parsing,
//! isotope physics, or ProFormA chemistry.
//!
//! This is deliberately **not** a second algorithm: every command maps 1:1 onto
//! a core function, and `tests/cli_contract.rs` asserts the CLI's `result` JSON
//! is byte-identical to a direct core call.

mod args;
mod commands;
mod output;
mod unimod_db;

use std::process::ExitCode;

use clap::Parser;

use args::Args;
use commands::DispatchOutput;

fn main() -> ExitCode {
    let args = match Args::try_parse() {
        Ok(args) => args,
        Err(err) => {
            use clap::error::ErrorKind;
            return match err.kind() {
                ErrorKind::DisplayHelp | ErrorKind::DisplayVersion => {
                    // `--help` is the interface; print it plainly.
                    let _ = err.print();
                    ExitCode::SUCCESS
                }
                _ => {
                    // Usage errors are a stable JSON shape on stderr, exit 2.
                    output::emit_usage_error(&err.to_string());
                    ExitCode::from(2)
                }
            };
        }
    };

    let started = std::time::Instant::now();
    match commands::dispatch(&args) {
        Ok(DispatchOutput::Json(mut envelope)) => {
            envelope.elapsed_ms = started.elapsed().as_millis();
            output::emit_ok(&envelope);
            ExitCode::SUCCESS
        }
        Ok(DispatchOutput::Mgf(text)) => {
            // A Mascot Generic Format spectrum list (scans --out-format mgf).
            print!("{text}");
            ExitCode::SUCCESS
        }
        Ok(DispatchOutput::Tab(text)) => {
            // A standard tab-separated document (fragment-theoretical
            // --out-format mztab): the mzTab-M molecular table is the whole
            // stdout so external consumers can read it directly.
            print!("{text}");
            ExitCode::SUCCESS
        }
        Err(error) => {
            output::emit_json_error(&error.to_json());
            ExitCode::FAILURE
        }
    }
}
