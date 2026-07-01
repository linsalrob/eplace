"""Header parsers for non-standard reference FASTA databases."""

import re
from typing import Optional, Tuple


def is_missing_taxid(value: Optional[str]) -> bool:
    """Return True when BLAST did not provide a usable taxid."""
    return value is None or value.strip() in {"", "0", "N/A", "NA", "-"}


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

    taxid_match = re.search(r"\|taxid=(\d+)", header)
    taxid = taxid_match.group(1) if taxid_match else ""

    desc_match = re.search(r"^" + re.escape(seq_id) + r"\s+([^\[]*)", header)
    description = desc_match.group(1).strip() if desc_match else ""

    return seq_id, description, taxid


def recover_taxid_from_subject_title(subject_title: str) -> Tuple[str, str]:
    """Recover subject ID and taxid from an NBDL BLAST stitle field."""
    seq_id, _description, taxid = parse_nbdl_header(subject_title)
    return seq_id, taxid
