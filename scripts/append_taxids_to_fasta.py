#!/usr/bin/env python3

"""
Append taxids to custom reference FASTA headers using organism/species names.

This utility is intended for custom reference databases such as CSIRO/NBDL,
BOLD, institutional FASTA exports, or mixed curated databases where sequence
headers contain organism names but may not already contain NCBI taxids.

The script tries to:

1. Extract a species/organism name from flexible FASTA header formats.
2. Map the exact species name to an NCBI taxid using names2id.tsv.
3. If the exact species is absent from NCBI, fall back to genus-level taxid
   when the genus exists in names2id.tsv.
4. Preserve the original custom organism name in the header using
   [eplace_original_organism=...] so downstream tools can still recover the
   informative custom species label.
5. Write a detailed TSV report recording exact, genus-fallback, and unmapped
   records.

Expected names2id.tsv format:

    scientific_name<TAB>taxid

Example:

    Chimaera fulva<TAB>765187
    Chimaera<TAB>30331

Important:

- The appended |taxid= value should always be a real taxid present in your
  taxonomy dump. Do not append fake custom IDs as |taxid= unless your taxonomy
  database has been extended to include them.
- For species absent from NCBI, the script appends the nearest valid genus
  taxid and preserves the original species name in header metadata.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


RANKED_NAME_PATTERNS = [
    # [organism=Chimaera fulva]
    ("organism_bracket", re.compile(r"\[organism=([^\]]+)\]", re.IGNORECASE)),

    # [species=Chimaera fulva] or [species=https://.../Chimaera_fulva]
    ("species_bracket", re.compile(r"\[species=([^\]]+)\]", re.IGNORECASE)),

    # organism=Chimaera fulva; or organism="Chimaera fulva"
    ("organism_key", re.compile(r"\borganism=['\"]?([^;,\]\[]+)['\"]?", re.IGNORECASE)),

    # species=Chimaera fulva; or species=https://.../Chimaera_fulva
    ("species_key", re.compile(r"\bspecies=['\"]?([^;,\]\[]+)['\"]?", re.IGNORECASE)),
]

BAD_THIRD_WORDS = {
    "mitochondrially",
    "mitochondrial",
    "encoded",
    "ribosomal",
    "rrna",
    "rna",
    "gene",
    "complete",
    "partial",
    "voucher",
    "isolate",
    "specimen",
    "mitogenome",
    "genome",
}

HEADER_METADATA_KEYS = [
    "eplace_original_organism",
    "eplace_taxid_assignment",
    "eplace_taxid_rank",
    "eplace_taxid_name",
]


@dataclass
class NameMatch:
    status: str
    candidate_name: str
    assigned_taxid: str
    assigned_rank: str
    assigned_name: str
    extraction_method: str
    custom_label: str = ""


def clean_candidate_name(value: str) -> str:
    """
    Clean a candidate organism/species string.

    Handles URL endings like https://.../Chimaera_fulva, underscores, quotes,
    and repeated whitespace.
    """
    value = value.strip().strip('"').strip("'")

    if "://" in value:
        value = value.rstrip("/").split("/")[-1]

    value = value.replace("_", " ")
    value = re.sub(r"\s+", " ", value).strip()
    return value


def normalize_name_for_lookup(name: str) -> str:
    """Normalize names for case-insensitive dictionary matching."""
    name = clean_candidate_name(name)
    name = re.sub(r"\s+", " ", name).strip()
    return name.lower()


def extract_binomial_or_trinomial(text: str) -> Optional[str]:
    """
    Extract a plausible Latin binomial or trinomial from free text.

    This is a fallback for headers without explicit [organism=...] metadata.
    """
    cleaned = re.sub(r"\[[^\]]+\]", " ", text)
    cleaned = cleaned.replace("_", " ")
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

    # Prefer trinomial when the third word does not look like gene text.
    tri = re.search(
        r"\b([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)\b",
        cleaned,
    )
    if tri:
        third = tri.group(3).lower()
        if third not in BAD_THIRD_WORDS:
            return f"{tri.group(1)} {tri.group(2)} {tri.group(3)}"

    bi = re.search(
        r"\b([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)\b",
        cleaned,
    )
    if bi:
        return f"{bi.group(1)} {bi.group(2)}"

    return None


def extract_name_from_header(header: str) -> tuple[Optional[str], str]:
    """Return candidate organism/species name and extraction method."""
    h = header.lstrip(">").strip()

    for method, pattern in RANKED_NAME_PATTERNS:
        match = pattern.search(h)
        if not match:
            continue

        candidate = clean_candidate_name(match.group(1))
        if not candidate:
            continue

        # For URL-like species fields this should recover a clean binomial.
        binomial = extract_binomial_or_trinomial(candidate) or candidate
        return binomial, method

    free_candidate = extract_binomial_or_trinomial(h)
    if free_candidate:
        return free_candidate, "free_text_binomial"

    return None, "unresolved"


def genus_from_name(name: str) -> Optional[str]:
    """Return genus-like first token from a scientific name."""
    parts = clean_candidate_name(name).split()
    if not parts:
        return None

    genus = parts[0]
    if re.match(r"^[A-Z][a-zA-Z-]+$", genus):
        return genus
    return None


def load_names2id(path: Path) -> dict[str, str]:
    """Load names2id.tsv as a case-insensitive name -> taxid dictionary."""
    lookup: dict[str, str] = {}

    with path.open() as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue

            name = row[0].strip()
            taxid = row[1].strip()

            if not name or not taxid:
                continue

            lookup[normalize_name_for_lookup(name)] = taxid

    return lookup


def existing_taxid(header: str) -> Optional[str]:
    """Return existing |taxid=123 value if present."""
    match = re.search(r"\|taxid=(\d+)", header)
    if match:
        return match.group(1)
    return None


def strip_existing_taxid(header: str) -> str:
    """Remove existing terminal |taxid=123 from a header."""
    return re.sub(r"\s*\|taxid=\d+\s*$", "", header.rstrip())


def strip_eplace_metadata(header: str) -> str:
    """Remove existing ePLACE taxid-assignment metadata fields."""
    cleaned = header.rstrip()
    for key in HEADER_METADATA_KEYS:
        cleaned = re.sub(rf"\s*\[{re.escape(key)}=[^\]]*\]", "", cleaned)
    return cleaned


def bracket_escape(value: str) -> str:
    """Make a value safe for bracket-style FASTA metadata."""
    return value.replace("]", ")").replace("[", "(").strip()


def append_assignment_metadata(header: str, match: NameMatch) -> str:
    """
    Append ePLACE metadata preserving the original custom organism name.

    This is deliberately placed in the description area, not the first sequence
    token, so BLAST -parse_seqids behaviour remains stable.
    """
    if not match.candidate_name:
        return header

    fields = [
        f"[eplace_original_organism={bracket_escape(match.candidate_name)}]",
        f"[eplace_taxid_assignment={match.status}]",
        f"[eplace_taxid_rank={match.assigned_rank or 'unmapped'}]",
    ]

    if match.assigned_name:
        fields.append(f"[eplace_taxid_name={bracket_escape(match.assigned_name)}]")

    if match.custom_label:
        fields.append(f"[eplace_custom_label={bracket_escape(match.custom_label)}]")

    return f"{header} {' '.join(fields)}"


def append_taxid_to_header(
    header: str,
    taxid: Optional[str],
    *,
    replace_existing: bool = False,
) -> str:
    """Append |taxid=<taxid> to a FASTA header."""
    header = header.rstrip()

    if existing_taxid(header) and not replace_existing:
        return header

    if existing_taxid(header) and replace_existing:
        header = strip_existing_taxid(header)

    if taxid:
        return f"{header} |taxid={taxid}"

    return header


def resolve_taxid(
    candidate_name: Optional[str],
    names2id: dict[str, str],
    *,
    custom_taxon_prefix: Optional[str] = None,
    custom_taxa_seen: Optional[dict[str, str]] = None,
) -> NameMatch:
    """
    Resolve candidate name to exact species taxid, genus taxid, or unmapped.

    The genus fallback is intentionally conservative: the original species name
    is preserved in ePLACE metadata while |taxid= points to a real NCBI genus.
    """
    if not candidate_name:
        return NameMatch(
            status="unmapped_no_candidate_name",
            candidate_name="",
            assigned_taxid="",
            assigned_rank="",
            assigned_name="",
            extraction_method="unresolved",
        )

    normalized_candidate = normalize_name_for_lookup(candidate_name)
    exact_taxid = names2id.get(normalized_candidate)
    if exact_taxid:
        return NameMatch(
            status="mapped_exact",
            candidate_name=candidate_name,
            assigned_taxid=exact_taxid,
            assigned_rank="species_or_name",
            assigned_name=candidate_name,
            extraction_method="",
        )

    genus = genus_from_name(candidate_name)
    if genus:
        genus_taxid = names2id.get(normalize_name_for_lookup(genus))
        if genus_taxid:
            return NameMatch(
                status="mapped_genus_fallback",
                candidate_name=candidate_name,
                assigned_taxid=genus_taxid,
                assigned_rank="genus",
                assigned_name=genus,
                extraction_method="",
            )

    custom_label = ""
    if custom_taxon_prefix and custom_taxa_seen is not None:
        key = normalized_candidate
        if key not in custom_taxa_seen:
            custom_taxa_seen[key] = f"{custom_taxon_prefix}{len(custom_taxa_seen) + 1}"
        custom_label = custom_taxa_seen[key]

    return NameMatch(
        status="unmapped",
        candidate_name=candidate_name,
        assigned_taxid="",
        assigned_rank="",
        assigned_name="",
        extraction_method="",
        custom_label=custom_label,
    )


def process_fasta(
    input_fasta: Path,
    output_fasta: Path,
    names2id: dict[str, str],
    report_tsv: Path,
    *,
    replace_existing: bool = False,
    refresh_eplace_metadata: bool = False,
    custom_taxon_prefix: Optional[str] = None,
) -> None:
    """Process FASTA headers and append usable taxids plus ePLACE metadata."""
    total = 0
    already_had_taxid = 0
    mapped_exact = 0
    mapped_genus = 0
    unresolved = 0
    custom_taxa_seen: dict[str, str] = {}

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)

    with input_fasta.open() as inp, output_fasta.open("w") as out, report_tsv.open("w") as rep:
        report_writer = csv.writer(rep, delimiter="\t")
        report_writer.writerow([
            "seq_id",
            "status",
            "candidate_name",
            "assigned_taxid",
            "assigned_rank",
            "assigned_name",
            "custom_label",
            "extraction_method",
            "original_header",
            "new_header",
        ])

        for line in inp:
            if not line.startswith(">"):
                out.write(line)
                continue

            total += 1
            original_header = line.rstrip("\n")
            seq_id = original_header.lstrip(">").split()[0]

            current_taxid = existing_taxid(original_header)
            candidate_name, method = extract_name_from_header(original_header)

            working_header = original_header
            if refresh_eplace_metadata:
                working_header = strip_eplace_metadata(working_header)

            if current_taxid and not replace_existing:
                already_had_taxid += 1

                match = NameMatch(
                    status="already_had_taxid",
                    candidate_name=candidate_name or "",
                    assigned_taxid=current_taxid,
                    assigned_rank="existing",
                    assigned_name=candidate_name or "",
                    extraction_method=method,
                )

                if refresh_eplace_metadata and candidate_name:
                    working_header = append_assignment_metadata(working_header, match)

                out.write(working_header + "\n")
                report_writer.writerow([
                    seq_id,
                    match.status,
                    match.candidate_name,
                    match.assigned_taxid,
                    match.assigned_rank,
                    match.assigned_name,
                    match.custom_label,
                    method,
                    original_header,
                    working_header,
                ])
                continue

            match = resolve_taxid(
                candidate_name,
                names2id,
                custom_taxon_prefix=custom_taxon_prefix,
                custom_taxa_seen=custom_taxa_seen,
            )
            match.extraction_method = method

            if match.status == "mapped_exact":
                mapped_exact += 1
            elif match.status == "mapped_genus_fallback":
                mapped_genus += 1
            else:
                unresolved += 1

            working_header = append_assignment_metadata(working_header, match)
            working_header = append_taxid_to_header(
                working_header,
                match.assigned_taxid or None,
                replace_existing=replace_existing,
            )

            out.write(working_header + "\n")
            report_writer.writerow([
                seq_id,
                match.status,
                match.candidate_name,
                match.assigned_taxid,
                match.assigned_rank,
                match.assigned_name,
                match.custom_label,
                method,
                original_header,
                working_header,
            ])

    print(f"Input records: {total}")
    print(f"Already had taxid: {already_had_taxid}")
    print(f"Mapped exact name: {mapped_exact}")
    print(f"Mapped by genus fallback: {mapped_genus}")
    print(f"Unresolved/no usable ancestor: {unresolved}")
    print(f"Wrote FASTA: {output_fasta}")
    print(f"Wrote report: {report_tsv}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Append NCBI taxids to custom FASTA headers using organism/species "
            "names, with genus-level fallback for taxa absent from NCBI."
        )
    )
    parser.add_argument(
        "-i",
        "--input-fasta",
        required=True,
        type=Path,
        help="Input FASTA file.",
    )
    parser.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        type=Path,
        help="Output FASTA file with |taxid= appended where possible.",
    )
    parser.add_argument(
        "-n",
        "--names2id",
        required=True,
        type=Path,
        help="Two-column TSV: scientific_name<TAB>taxid.",
    )
    parser.add_argument(
        "-r",
        "--report-tsv",
        required=True,
        type=Path,
        help="Output TSV report of exact, genus-fallback, and unmapped records.",
    )
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help="Replace existing |taxid= values instead of preserving them.",
    )
    parser.add_argument(
        "--refresh-eplace-metadata",
        action="store_true",
        help=(
            "Remove old [eplace_*] metadata fields and rewrite them. Useful when "
            "rerunning this script on a previously annotated FASTA."
        ),
    )
    parser.add_argument(
        "--custom-taxon-prefix",
        default="CUSTOM_",
        help=(
            "Prefix for unresolved custom taxa in the report/header metadata. "
            "These labels are not appended as |taxid=. Default: CUSTOM_"
        ),
    )

    args = parser.parse_args()

    names2id = load_names2id(args.names2id)

    process_fasta(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        names2id=names2id,
        report_tsv=args.report_tsv,
        replace_existing=args.replace_existing,
        refresh_eplace_metadata=args.refresh_eplace_metadata,
        custom_taxon_prefix=args.custom_taxon_prefix,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
