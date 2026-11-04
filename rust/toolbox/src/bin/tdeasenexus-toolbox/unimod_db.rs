//! SQLite reader for the bundled UniMod database (D12-R). This is the only
//! place the CLI touches SQL; the search/composition/solve logic lives in the
//! lib (`unimod.rs`) and operates on the records loaded here.
//!
//! The schema is owned by `scripts/build-unimod-database.mjs` (schema v2).
//! Open read-only so a bug can never corrupt the shipped asset.

use rusqlite::{Connection, OpenFlags};
use sha2::{Digest, Sha256};
use tdeasenexus_toolbox::{
    NameKind, UniModAminoAcid, UniModNameEntry, UniModNeutralLoss, UniModRecord, UniModSourceMeta,
    UniModSpecificity, USER_OFFSET,
};

const EMBEDDED_UNIMOD_DB: &[u8] = include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../system/reference/unimod.sqlite"
));

/// Resolve the reference database path. `--db` is a development/test override;
/// production uses the database compiled into this binary.
pub fn resolve_db_path(explicit: Option<&str>) -> Result<String, crate::commands::CommandError> {
    if let Some(path) = explicit {
        if !path.trim().is_empty() {
            return Ok(path.to_owned());
        }
    }
    materialize_embedded_db()
}

fn materialize_embedded_db() -> Result<String, crate::commands::CommandError> {
    let digest = format!("{:x}", Sha256::digest(EMBEDDED_UNIMOD_DB));
    let root = std::env::temp_dir().join("tdeasenexus-toolbox");
    std::fs::create_dir_all(&root).map_err(|error| crate::commands::CommandError::FileRead {
        path: root.to_string_lossy().into_owned(),
        error: error.to_string(),
    })?;
    let target = root.join(format!("unimod-{digest}.sqlite"));
    if !target.is_file() {
        let pending = root.join(format!("unimod-{digest}-{}.tmp", std::process::id()));
        match std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&pending)
        {
            Ok(mut file) => {
                use std::io::Write;
                file.write_all(EMBEDDED_UNIMOD_DB).map_err(|error| {
                    crate::commands::CommandError::FileRead {
                        path: pending.to_string_lossy().into_owned(),
                        error: error.to_string(),
                    }
                })?;
                file.sync_all()
                    .map_err(|error| crate::commands::CommandError::FileRead {
                        path: pending.to_string_lossy().into_owned(),
                        error: error.to_string(),
                    })?;
                if let Err(error) = std::fs::rename(&pending, &target) {
                    if !target.is_file() {
                        return Err(crate::commands::CommandError::FileRead {
                            path: target.to_string_lossy().into_owned(),
                            error: error.to_string(),
                        });
                    }
                    let _ = std::fs::remove_file(&pending);
                }
            }
            Err(error) if error.kind() == std::io::ErrorKind::AlreadyExists => {}
            Err(error) => {
                return Err(crate::commands::CommandError::FileRead {
                    path: pending.to_string_lossy().into_owned(),
                    error: error.to_string(),
                })
            }
        }
    }
    Ok(target.to_string_lossy().into_owned())
}

pub struct UnimodDb {
    conn: Connection,
}

impl UnimodDb {
    pub fn open(path: &str) -> Result<Self, crate::commands::CommandError> {
        let conn = Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY).map_err(
            |error| crate::commands::CommandError::FileRead {
                path: path.to_owned(),
                error: error.to_string(),
            },
        )?;
        Ok(Self { conn })
    }

    pub fn meta(&self) -> Result<UniModSourceMeta, crate::commands::CommandError> {
        let get = |key: &str| -> Result<String, crate::commands::CommandError> {
            self.conn
                .query_row(
                    "SELECT value FROM reference_metadata WHERE key = ?1",
                    [key],
                    |row| row.get::<_, String>(0),
                )
                .map_err(|error| crate::commands::CommandError::FileRead {
                    path: "unimod.sqlite metadata".to_owned(),
                    error: error.to_string(),
                })
        };
        Ok(UniModSourceMeta {
            source_sha256: get("source_sha256")?,
            schema_version: get("schema_version")?.parse().unwrap_or_default(),
            release_major: get("unimod_release_major")?,
            release_minor: get("unimod_release_minor")?,
        })
    }

    pub fn all_records(&self) -> Result<Vec<UniModRecord>, crate::commands::CommandError> {
        let mut stmt = self
            .conn
            .prepare("SELECT record_id, title, full_name, approved, delta_mono, delta_average, composition FROM modifications ORDER BY record_id")
            .map_err(db_error)?;
        let rows = stmt
            .query_map([], |row| {
                Ok(UniModRecord {
                    record_id: row.get(0)?,
                    title: row.get(1)?,
                    full_name: row.get(2)?,
                    approved: row.get::<_, i64>(3)? != 0,
                    delta_mono: row.get(4)?,
                    delta_average: row.get(5)?,
                    composition: row.get(6)?,
                    source: "unimod".to_owned(),
                })
            })
            .map_err(db_error)?;
        rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
    }

    pub fn names(&self) -> Result<Vec<UniModNameEntry>, crate::commands::CommandError> {
        let mut stmt = self
            .conn
            .prepare("SELECT modification_id, name, kind FROM modification_names")
            .map_err(db_error)?;
        let rows = stmt
            .query_map([], |row| {
                let kind: String = row.get(2)?;
                Ok(UniModNameEntry {
                    modification_id: row.get(0)?,
                    name: row.get(1)?,
                    kind: match kind.as_str() {
                        "title" => NameKind::Title,
                        "full_name" => NameKind::FullName,
                        _ => NameKind::AltName,
                    },
                })
            })
            .map_err(db_error)?;
        rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
    }

    pub fn all_specificities(
        &self,
    ) -> Result<Vec<UniModSpecificity>, crate::commands::CommandError> {
        let mut stmt = self
            .conn
            .prepare("SELECT id, modification_id, site, position, classification, spec_group_id, hidden FROM specificities")
            .map_err(db_error)?;
        let rows = stmt
            .query_map([], |row| {
                Ok(UniModSpecificity {
                    id: row.get(0)?,
                    modification_id: row.get(1)?,
                    site: row.get(2)?,
                    position: row.get(3)?,
                    classification: row.get(4)?,
                    spec_group_id: row.get(5)?,
                    hidden: row.get::<_, i64>(6)? != 0,
                })
            })
            .map_err(db_error)?;
        rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
    }

    pub fn all_neutral_losses(
        &self,
    ) -> Result<Vec<UniModNeutralLoss>, crate::commands::CommandError> {
        let mut stmt = self
            .conn
            .prepare(
                "SELECT specificity_id, mono_mass, average_mass, composition FROM neutral_losses",
            )
            .map_err(db_error)?;
        let rows = stmt
            .query_map([], |row| {
                Ok(UniModNeutralLoss {
                    specificity_id: row.get(0)?,
                    mono_mass: row.get(1)?,
                    average_mass: row.get(2)?,
                    composition: row.get(3)?,
                })
            })
            .map_err(db_error)?;
        rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
    }

    pub fn amino_acids(&self) -> Result<Vec<UniModAminoAcid>, crate::commands::CommandError> {
        let mut stmt = self
            .conn
            .prepare(
                "SELECT title, three_letter, full_name, mono_mass, average_mass FROM amino_acids",
            )
            .map_err(db_error)?;
        let rows = stmt
            .query_map([], |row| {
                Ok(UniModAminoAcid {
                    title: row.get(0)?,
                    three_letter: row.get(1)?,
                    full_name: row.get(2)?,
                    mono_mass: row.get(3)?,
                    average_mass: row.get(4)?,
                })
            })
            .map_err(db_error)?;
        rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
    }
}

fn db_error(error: rusqlite::Error) -> crate::commands::CommandError {
    crate::commands::CommandError::FileRead {
        path: "unimod.sqlite".to_owned(),
        error: error.to_string(),
    }
}

// ── User modification store (D17) ──────────────────────────────────────────
//
// `user-mods.sqlite` is a writable User System store (`USER:n`) with the same
// projected schema as the official database. Its record ids are offset above
// the official space so merged records never collide across sources, and every
// record is tagged `source="user"` / `custom=true` downstream.

/// Load the user store's records (empty when the store does not exist yet).
pub fn load_user_records(path: &str) -> Result<Vec<UniModRecord>, crate::commands::CommandError> {
    if !std::path::Path::new(path).exists() {
        return Ok(Vec::new());
    }
    let conn =
        Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY).map_err(db_error)?;
    let mut stmt = conn
        .prepare("SELECT record_id, title, full_name, approved, delta_mono, delta_average, composition FROM modifications ORDER BY record_id")
        .map_err(db_error)?;
    let rows = stmt
        .query_map([], |row| {
            Ok(UniModRecord {
                record_id: row.get::<_, u32>(0)?,
                title: row.get(1)?,
                full_name: row.get(2)?,
                approved: row.get::<_, i64>(3)? != 0,
                delta_mono: row.get(4)?,
                delta_average: row.get(5)?,
                composition: row.get(6)?,
                source: "user".to_owned(),
            })
        })
        .map_err(db_error)?;
    rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
}

/// Load the user store's names (empty when the store does not exist yet).
pub fn load_user_names(path: &str) -> Result<Vec<UniModNameEntry>, crate::commands::CommandError> {
    if !std::path::Path::new(path).exists() {
        return Ok(Vec::new());
    }
    let conn =
        Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY).map_err(db_error)?;
    let mut stmt = conn
        .prepare("SELECT modification_id, name, kind FROM modification_names")
        .map_err(db_error)?;
    let rows = stmt
        .query_map([], |row| {
            let kind: String = row.get(2)?;
            Ok(UniModNameEntry {
                modification_id: row.get(0)?,
                name: row.get(1)?,
                kind: match kind.as_str() {
                    "title" => NameKind::Title,
                    "full_name" => NameKind::FullName,
                    _ => NameKind::AltName,
                },
            })
        })
        .map_err(db_error)?;
    rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
}

/// Load the user store's specificities (empty when the store does not exist yet).
pub fn load_user_specificities(
    path: &str,
) -> Result<Vec<UniModSpecificity>, crate::commands::CommandError> {
    if !std::path::Path::new(path).exists() {
        return Ok(Vec::new());
    }
    let conn =
        Connection::open_with_flags(path, OpenFlags::SQLITE_OPEN_READ_ONLY).map_err(db_error)?;
    let mut stmt = conn
        .prepare("SELECT id, modification_id, site, position, classification, spec_group_id, hidden FROM specificities")
        .map_err(db_error)?;
    let rows = stmt
        .query_map([], |row| {
            Ok(UniModSpecificity {
                id: row.get(0)?,
                modification_id: row.get(1)?,
                site: row.get(2)?,
                position: row.get(3)?,
                classification: row.get(4)?,
                spec_group_id: row.get(5)?,
                hidden: row.get::<_, i64>(6)? != 0,
            })
        })
        .map_err(db_error)?;
    rows.collect::<Result<Vec<_>, _>>().map_err(db_error)
}

/// Add one user modification (`USER:n`). Returns the new user id (n).
pub fn add_user_modification(
    path: &str,
    title: &str,
    full_name: &str,
    formula: &str,
    mono_mass: f64,
    average_mass: f64,
    sites: &[String],
) -> Result<u32, crate::commands::CommandError> {
    let conn = Connection::open(path).map_err(db_error)?;
    conn.execute_batch(
        "CREATE TABLE IF NOT EXISTS modifications (
            record_id INTEGER PRIMARY KEY,
            title TEXT NOT NULL,
            full_name TEXT NOT NULL,
            approved INTEGER NOT NULL,
            delta_mono REAL NOT NULL,
            delta_average REAL NOT NULL,
            composition TEXT NOT NULL DEFAULT '',
            username_of_poster TEXT NOT NULL DEFAULT '',
            group_of_poster TEXT NOT NULL DEFAULT '',
            date_time_posted TEXT NOT NULL DEFAULT '',
            date_time_modified TEXT NOT NULL DEFAULT ''
        );
        CREATE TABLE IF NOT EXISTS modification_names (
            modification_id INTEGER NOT NULL,
            name TEXT NOT NULL,
            kind TEXT NOT NULL CHECK (kind IN ('title','full_name','alt_name'))
        );
        CREATE TABLE IF NOT EXISTS specificities (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            modification_id INTEGER NOT NULL,
            site TEXT NOT NULL,
            position TEXT NOT NULL,
            classification TEXT NOT NULL,
            spec_group_id INTEGER,
            hidden INTEGER NOT NULL
        );
        CREATE TABLE IF NOT EXISTS reference_metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL);",
    )
    .map_err(db_error)?;

    // The next user id is the highest existing user id + 1.
    let next_user: u32 = conn
        .query_row(
            "SELECT COALESCE(MAX(record_id) - ?1, 0) + 1 FROM modifications WHERE record_id >= ?1",
            [USER_OFFSET],
            |row| row.get(0),
        )
        .map_err(db_error)?;
    let record_id = USER_OFFSET + next_user;

    conn.execute(
        "INSERT INTO modifications (record_id, title, full_name, approved, delta_mono, delta_average, composition) VALUES (?1,?2,?3,1,?4,?5,?6)",
        rusqlite::params![record_id, title, full_name, mono_mass, average_mass, formula],
    )
    .map_err(db_error)?;
    conn.execute(
        "INSERT INTO modification_names (modification_id, name, kind) VALUES (?1,?2,'title')",
        rusqlite::params![record_id, title],
    )
    .map_err(db_error)?;
    for site in sites {
        conn.execute(
            "INSERT INTO specificities (modification_id, site, position, classification, hidden) VALUES (?1,?2,'','Anywhere',0)",
            rusqlite::params![record_id, site],
        )
        .map_err(db_error)?;
    }
    Ok(next_user)
}
