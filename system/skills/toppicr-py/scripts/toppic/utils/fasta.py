"""
FASTA file utilities: reading, parsing, and classifying protein accessions.

Uses BioPython for FASTA I/O. Provides helpers for the TopPIC pipeline
to identify decoy entries and classify annotation types.
"""

import gzip
import re
from pathlib import Path
from typing import Optional


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file and return {accession: sequence} dictionary.

    Supports both plain-text ``.fasta`` and gzip-compressed ``.fasta.gz``
    files. Standard-library implementation (no BioPython), so the skill runs
    from a minimal wheelhouse.

    Parameters
    ----------
    path : str | Path
        Path to .fasta or .fasta.gz file.

    Returns
    -------
    dict[str, str]
        Keys are the full accession header (first word after ``>``).
        Values are the protein sequences as strings.
    """
    records = {}
    path = Path(path)

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(str(path), "rt") as fh:
        accession = None
        seq_parts = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = _clean_sequence("".join(seq_parts))
                accession = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if accession is not None:
            records[accession] = _clean_sequence("".join(seq_parts))
    return records


def _clean_sequence(seq: str) -> str:
    """Replace ambiguous amino acids as done in the R code."""
    return seq.replace("X", "A").replace("B", "D").replace("Z", "E").replace("J", "I")


def is_decoy(accession: str) -> bool:
    """Check whether a protein accession is a TopPIC decoy entry."""
    return accession.startswith("DECOY") or accession.startswith("DECOY_")


def classify_annotation(accession: str) -> Optional[str]:
    """Classify a UniProt accession into annotation type.

    Parameters
    ----------
    accession : str
        Protein accession header (the full first word, e.g.
        ``sp|P17096|HMGA1_HUMAN`` or ``tr|A0A0A0MS03|...``).

    Returns
    -------
    str or None
        ``"SwissProt"``, ``"TrEMBL"``, ``"VarSplic"``, or ``None``.
    """
    decoy = re.escape("DECOY_")
    pattern = r"^(DECOY_)?sp\|([^-]*)\|.*"
    if re.match(rf"^(?:{decoy})?sp\|[^-]+-\d+\|.*", accession):
        return "VarSplic"
    if re.match(rf"^(?:{decoy})?tr.*", accession):
        return "TrEMBL"
    if re.match(rf"^(?:{decoy})?sp\|[^-]*\|.*", accession):
        return "SwissProt"
    return None


def extract_uniprot_acc(accession_header: str) -> str:
    """Extract UniProt accession from a full accession header.

    The header format is typically ``sp|P17096|HMGA1_HUMAN`` or
    ``tr|A0A0A0MS03|...`` — the second pipe-delimited field is the
    actual accession.

    Parameters
    ----------
    accession_header : str
        Full protein accession header line.

    Returns
    -------
    str
        The UniProt accession string (e.g. ``"P17096"``).
    """
    parts = accession_header.split("|")
    if len(parts) >= 2:
        return parts[1]
    return accession_header


def extract_gene(protein_description: str) -> str:
    """Extract the gene name from a protein description string.

    Looks for the ``GN=`` tag in the description (e.g.
    ``... GN=HMGA1 PE=1 ...``).

    Parameters
    ----------
    protein_description : str
        The protein description field from TopPIC output.

    Returns
    -------
    str
        The gene name if found, else ``""``.
    """
    m = re.search(r"GN=(\S+)", protein_description)
    return m.group(1) if m else ""
