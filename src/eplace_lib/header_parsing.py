"""Header parsers for non-standard reference FASTA databases."""

import re
from typing import Optional, Tuple


MISSING_TAXID_VALUES = {"", "0", "N/A", "NA", "-", "NONE", "NULL"}


def is_missing_taxid(value: Optional[str]) -> bool:
    """Return True when BLAST did not provide a usable taxid."""
    if value is None:
        return True
    return value.strip().upper() in MISSING_TAXID_VALUES


def extract_taxid_from_text(text: str) -> str:
    """Extract taxid=123 from arbitrary header or BLAST title text."""
    if not text:
        return ""
    match = re.search(r"(?:\||\[|\s)taxid=(\d+)(?:\]|\s|$)", text, re.IGNORECASE)
    return match.group(1) if match else ""


def parse_nbdl_header(header: str) -> Tuple[str, str, str]:
    """Parse an NBDL-style FASTA header.

    Returns:
        sequence_id, description, taxid
    """
    if not header:
        return "", "", ""

    header = header.lstrip(">")

    seq_id_match = re.match(r"([^\s\[]+)", header)
    seq_id = seq_id_match.group(1) if seq_id_match else header

    taxid = extract_taxid_from_text(header)

    desc_match = re.search(r"^" + re.escape(seq_id) + r"\s+([^\[]*)", header)
    description = desc_match.group(1).strip() if desc_match else ""

    return seq_id, description, taxid


def recover_taxid_from_subject_title(subject_title: str) -> Tuple[str, str]:
    """Recover subject ID and taxid from an NBDL BLAST stitle field."""
    seq_id, _description, taxid = parse_nbdl_header(subject_title)
    if not taxid:
        taxid = extract_taxid_from_text(subject_title)
    return seq_id, taxid
