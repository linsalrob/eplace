#!/usr/bin/env python3

"""
Append usable NCBI taxids to custom reference FASTA headers while preserving
original sequence identifiers.

This utility is intended for custom reference databases such as CSIRO/NBDL,
BOLD, institutional FASTA exports, or mixed curated databases where sequence
headers contain organism names but may not already contain NCBI taxids.

The FASTA output deliberately preserves the original first FASTA token. For an
NBDL record, the output header remains NBDL-native and only receives a terminal
|taxid=<taxid> annotation when a usable NCBI taxid or conservative NCBI ancestor
can be found.

Example output header:

    >NBDL-HK3QATBVPPF8GP.v1.mt|15366-17013|+|MT-RNR2 [species=...] ... |taxid=765187

Resolution ladder:

1. Extract organism/species name from flexible FASTA header formats.
2. Map the exact name to a taxid using names2id.tsv.
3. If exact name is absent, try binomial fallback for trinomial names.
4. If exact/binomial name is absent, fall back to genus-level taxid when genus exists.
5. If genus is absent and --taxonomy-hierarchy is supplied, use that hierarchy
   to find family/order/class/superclass for the species or genus, then map the
   first available ancestor name to names2id.tsv.
6. If no valid NCBI ancestor can be found, preserve the custom organism name in
   the mapping tables, assign a stable custom label, and do not append fake
   |taxid= values.

The script writes:

1. Rewritten FASTA with original headers preserved and usable |taxid= annotations.
2. Report TSV describing exact, fallback, and unresolved mappings.
3. Reference mapping TSV intended for downstream ePLACE use. This file preserves
   the NBDL sequence ID/header, extracted organism, resolved NCBI taxid if any,
   fallback rank/name, hierarchy source, and custom unresolved label.

Expected names2id.tsv format:

    scientific_name<TAB>taxid

Expected --taxonomy-hierarchy CSV columns, case-insensitive:

    Species, Genus, Subfamily, Family, Order, Class, SuperClass

Important:

- The appended |taxid= value should always be a real taxid present in your
  taxonomy dump. Do not append fake custom IDs as |taxid= unless your taxonomy
  database has been extended to include them.
- For custom taxa absent from NCBI, the script preserves the original name in the
  mapping files even when |taxid= points to a conservative genus/family/order/
  class/superclass ancestor.
- This script does not sanitize or replace NBDL sequence IDs by default because
  the current ePLACE NBDL workflow recovers taxids from BLAST stitle/header text.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


RANKED_NAME_PATTERNS = [
    ("organism_bracket", re.compile(r"\[organism=([^\]]+)\]", re.IGNORECASE)),
    ("species_bracket", re.compile(r"\[species=([^\]]+)\]", re.IGNORECASE)),
    ("organism_key", re.compile(r"\borganism=['\"]?([^;,\]\[]+)['\"]?", re.IGNORECASE)),
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

DEFAULT_HIERARCHY_FALLBACK_RANKS = ["Family", "Order", "Class", "SuperClass"]


@dataclass
class HierarchyRecord:
    species: str = ""
    genus: str = ""
    subfamily: str = ""
    family: str = ""
    order: str = ""
    class_name: str = ""
    superclass: str = ""

    def value_for_rank(self, rank: str) -> str:
        rank_key = normalize_column_name(rank)
        if rank_key == "species":
            return self.species
        if rank_key == "genus":
            return self.genus
        if rank_key == "subfamily":
            return self.subfamily
        if rank_key == "family":
            return self.family
        if rank_key == "order":
            return self.order
        if rank_key == "class":
            return self.class_name
        if rank_key == "superclass":
            return self.superclass
        return ""


@dataclass
class TaxonomyHierarchy:
    by_species: dict[str, HierarchyRecord] = field(default_factory=dict)
    by_genus: dict[str, HierarchyRecord] = field(default_factory=dict)

    def lookup(self, candidate_name: str) -> Optional[HierarchyRecord]:
        normalized = normalize_name_for_lookup(candidate_name)
        if normalized in self.by_species:
            return self.by_species[normalized]

        genus = genus_from_name(candidate_name)
        if genus:
            return self.by_genus.get(normalize_name_for_lookup(genus))

        return None


@dataclass
class NameMatch:
    status: str
    candidate_name: str
    assigned_taxid: str
    assigned_rank: str
    assigned_name: str
    extraction_method: str
    hierarchy_match: str = ""
    hierarchy_family: str = ""
    hierarchy_order: str = ""
    hierarchy_class: str = ""
    hierarchy_superclass: str = ""
    custom_label: str = ""

    @property
    def has_usable_taxid(self) -> bool:
        return bool(self.assigned_taxid and self.assigned_taxid != "0")


def normalize_column_name(value: str) -> str:
    """Normalize CSV column names for case-insensitive matching."""
    return re.sub(r"[^a-z0-9]", "", value.lower())


def clean_candidate_name(value: str) -> str:
    """Clean a candidate organism/species string."""
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
    """Extract a plausible Latin binomial or trinomial from free text."""
    cleaned = re.sub(r"\[[^\]]+\]", " ", text)
    cleaned = cleaned.replace("_", " ")
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

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


def binomial_from_name(name: str) -> Optional[str]:
    """Return Genus species from a longer scientific name when possible."""
    parts = clean_candidate_name(name).split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
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


def sequence_id_from_header(header: str) -> str:
    """Return the original first whitespace-delimited FASTA token without leading >."""
    return header.lstrip(">").strip().split()[0]


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


def row_value(row: dict[str, str], column_map: dict[str, str], desired_column: str) -> str:
    """Return a CSV row value by case-insensitive normalized column name."""
    actual = column_map.get(normalize_column_name(desired_column))
    if not actual:
        return ""
    return row.get(actual, "").strip()


def load_taxonomy_hierarchy(path: Optional[Path]) -> TaxonomyHierarchy:
    """Load an optional fish taxonomy hierarchy CSV."""
    hierarchy = TaxonomyHierarchy()

    if path is None:
        return hierarchy

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return hierarchy

        column_map = {normalize_column_name(name): name for name in reader.fieldnames}

        for row in reader:
            record = HierarchyRecord(
                species=clean_candidate_name(row_value(row, column_map, "Species")),
                genus=clean_candidate_name(row_value(row, column_map, "Genus")),
                subfamily=clean_candidate_name(row_value(row, column_map, "Subfamily")),
                family=clean_candidate_name(row_value(row, column_map, "Family")),
                order=clean_candidate_name(row_value(row, column_map, "Order")),
                class_name=clean_candidate_name(row_value(row, column_map, "Class")),
                superclass=clean_candidate_name(row_value(row, column_map, "SuperClass")),
            )

            if record.species:
                hierarchy.by_species[normalize_name_for_lookup(record.species)] = record

            if record.genus:
                genus_key = normalize_name_for_lookup(record.genus)
                hierarchy.by_genus.setdefault(genus_key, record)

    return hierarchy


def existing_taxid(header: str) -> Optional[str]:
    """Return existing |taxid=123 value if present."""
    match = re.search(r"\|taxid=(\d+)", header)
    if match:
        return match.group(1)
    return None


def strip_existing_taxid(header: str) -> str:
    """Remove existing terminal |taxid=123 from a header."""
    return re.sub(r"\s*\|taxid=\d+\s*$", "", header.rstrip())


def append_taxid_to_header(header: str, taxid: Optional[str], *, replace_existing: bool = False) -> str:
    """Append |taxid=<taxid> to the end of a FASTA header without changing anything else."""
    header = header.rstrip()

    if existing_taxid(header) and not replace_existing:
        return header

    if existing_taxid(header) and replace_existing:
        header = strip_existing_taxid(header)

    if taxid:
        return f"{header} |taxid={taxid}"

    return header


def hierarchy_details(record: Optional[HierarchyRecord]) -> dict[str, str]:
    if record is None:
        return {
            "hierarchy_match": "",
            "hierarchy_family": "",
            "hierarchy_order": "",
            "hierarchy_class": "",
            "hierarchy_superclass": "",
        }

    return {
        "hierarchy_match": record.species or record.genus,
        "hierarchy_family": record.family,
        "hierarchy_order": record.order,
        "hierarchy_class": record.class_name,
        "hierarchy_superclass": record.superclass,
    }


def try_hierarchy_fallback(
    candidate_name: str,
    names2id: dict[str, str],
    hierarchy: TaxonomyHierarchy,
    fallback_ranks: list[str],
) -> Optional[NameMatch]:
    """Try mapping candidate species/genus to an ancestor using taxonomy hierarchy."""
    record = hierarchy.lookup(candidate_name)
    if record is None:
        return None

    details = hierarchy_details(record)

    for rank in fallback_ranks:
        ancestor_name = record.value_for_rank(rank)
        if not ancestor_name:
            continue

        ancestor_taxid = names2id.get(normalize_name_for_lookup(ancestor_name))
        if not ancestor_taxid:
            continue

        normalized_rank = normalize_column_name(rank)
        status = f"mapped_{normalized_rank}_hierarchy_fallback"

        return NameMatch(
            status=status,
            candidate_name=candidate_name,
            assigned_taxid=ancestor_taxid,
            assigned_rank=rank.lower(),
            assigned_name=ancestor_name,
            extraction_method="",
            hierarchy_match=details["hierarchy_match"],
            hierarchy_family=details["hierarchy_family"],
            hierarchy_order=details["hierarchy_order"],
            hierarchy_class=details["hierarchy_class"],
            hierarchy_superclass=details["hierarchy_superclass"],
        )

    # Hierarchy knew the organism/genus, but none of its chosen ancestors were in names2id.
    return NameMatch(
        status="unmapped_hierarchy_match_no_ncbi_ancestor",
        candidate_name=candidate_name,
        assigned_taxid="",
        assigned_rank="",
        assigned_name="",
        extraction_method="",
        hierarchy_match=details["hierarchy_match"],
        hierarchy_family=details["hierarchy_family"],
        hierarchy_order=details["hierarchy_order"],
        hierarchy_class=details["hierarchy_class"],
        hierarchy_superclass=details["hierarchy_superclass"],
    )


def assign_custom_label(
    candidate_name: str,
    custom_taxon_prefix: Optional[str],
    custom_taxa_seen: Optional[dict[str, str]],
) -> str:
    if not custom_taxon_prefix or custom_taxa_seen is None:
        return ""

    key = normalize_name_for_lookup(candidate_name)
    if key not in custom_taxa_seen:
        custom_taxa_seen[key] = f"{custom_taxon_prefix}{len(custom_taxa_seen) + 1}"
    return custom_taxa_seen[key]


def resolve_taxid(
    candidate_name: Optional[str],
    names2id: dict[str, str],
    *,
    taxonomy_hierarchy: Optional[TaxonomyHierarchy] = None,
    hierarchy_fallback_ranks: Optional[list[str]] = None,
    custom_taxon_prefix: Optional[str] = None,
    custom_taxa_seen: Optional[dict[str, str]] = None,
) -> NameMatch:
    """Resolve candidate name to exact, binomial, genus, hierarchy-derived ancestor, or unmapped."""
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

    binomial = binomial_from_name(candidate_name)
    if binomial and normalize_name_for_lookup(binomial) != normalized_candidate:
        binomial_taxid = names2id.get(normalize_name_for_lookup(binomial))
        if binomial_taxid:
            return NameMatch(
                status="mapped_binomial_fallback",
                candidate_name=candidate_name,
                assigned_taxid=binomial_taxid,
                assigned_rank="species",
                assigned_name=binomial,
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

    if taxonomy_hierarchy is not None:
        hierarchy_match = try_hierarchy_fallback(
            candidate_name=candidate_name,
            names2id=names2id,
            hierarchy=taxonomy_hierarchy,
            fallback_ranks=hierarchy_fallback_ranks or DEFAULT_HIERARCHY_FALLBACK_RANKS,
        )
        if hierarchy_match:
            if not hierarchy_match.has_usable_taxid:
                hierarchy_match.custom_label = assign_custom_label(
                    candidate_name,
                    custom_taxon_prefix,
                    custom_taxa_seen,
                )
            return hierarchy_match

    return NameMatch(
        status="unmapped_no_ncbi_match",
        candidate_name=candidate_name,
        assigned_taxid="",
        assigned_rank="",
        assigned_name="",
        extraction_method="",
        custom_label=assign_custom_label(candidate_name, custom_taxon_prefix, custom_taxa_seen),
    )


def write_report_row(writer: csv.writer, values: dict[str, str]) -> None:
    writer.writerow([
        values["seq_id"],
        values["status"],
        values["candidate_name"],
        values["assigned_taxid"],
        values["assigned_rank"],
        values["assigned_name"],
        values["hierarchy_match"],
        values["hierarchy_family"],
        values["hierarchy_order"],
        values["hierarchy_class"],
        values["hierarchy_superclass"],
        values["custom_label"],
        values["extraction_method"],
        values["original_header"],
        values["new_header"],
    ])


def write_reference_map_row(writer: csv.writer, values: dict[str, str]) -> None:
    writer.writerow([
        values["seq_id"],
        values["candidate_name"],
        values["status"],
        values["assigned_taxid"],
        values["assigned_rank"],
        values["assigned_name"],
        values["hierarchy_match"],
        values["hierarchy_family"],
        values["hierarchy_order"],
        values["hierarchy_class"],
        values["hierarchy_superclass"],
        values["custom_label"],
        values["extraction_method"],
        values["original_header"],
        values["new_header"],
    ])


def process_fasta(
    input_fasta: Path,
    output_fasta: Path,
    names2id: dict[str, str],
    report_tsv: Path,
    reference_map_tsv: Path,
    *,
    taxonomy_hierarchy: Optional[TaxonomyHierarchy] = None,
    hierarchy_fallback_ranks: Optional[list[str]] = None,
    replace_existing: bool = False,
    custom_taxon_prefix: Optional[str] = None,
) -> None:
    """Process FASTA headers and append usable taxids only at the end of headers."""
    total = 0
    already_had_taxid = 0
    mapped_exact = 0
    mapped_binomial = 0
    mapped_genus = 0
    mapped_hierarchy = 0
    unresolved = 0
    custom_taxa_seen: dict[str, str] = {}

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)
    reference_map_tsv.parent.mkdir(parents=True, exist_ok=True)

    report_header = [
        "seq_id",
        "status",
        "candidate_name",
        "assigned_taxid",
        "assigned_rank",
        "assigned_name",
        "hierarchy_match",
        "hierarchy_family",
        "hierarchy_order",
        "hierarchy_class",
        "hierarchy_superclass",
        "custom_label",
        "extraction_method",
        "original_header",
        "new_header",
    ]

    with (
        input_fasta.open() as inp,
        output_fasta.open("w") as out,
        report_tsv.open("w") as rep,
        reference_map_tsv.open("w") as refmap,
    ):
        report_writer = csv.writer(rep, delimiter="\t")
        refmap_writer = csv.writer(refmap, delimiter="\t")
        report_writer.writerow(report_header)
        refmap_writer.writerow(report_header)

        for line in inp:
            if not line.startswith(">"):
                out.write(line)
                continue

            total += 1
            original_header = line.rstrip("\n")
            seq_id = sequence_id_from_header(original_header)
            current_taxid = existing_taxid(original_header)
            candidate_name, method = extract_name_from_header(original_header)

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
                new_header = original_header
            else:
                match = resolve_taxid(
                    candidate_name,
                    names2id,
                    taxonomy_hierarchy=taxonomy_hierarchy,
                    hierarchy_fallback_ranks=hierarchy_fallback_ranks,
                    custom_taxon_prefix=custom_taxon_prefix,
                    custom_taxa_seen=custom_taxa_seen,
                )
                match.extraction_method = method

                if match.status == "mapped_exact":
                    mapped_exact += 1
                elif match.status == "mapped_binomial_fallback":
                    mapped_binomial += 1
                elif match.status == "mapped_genus_fallback":
                    mapped_genus += 1
                elif "hierarchy_fallback" in match.status:
                    mapped_hierarchy += 1
                else:
                    unresolved += 1

                new_header = append_taxid_to_header(
                    original_header,
                    match.assigned_taxid or None,
                    replace_existing=replace_existing,
                )

            out.write(new_header + "\n")

            values = {
                "seq_id": seq_id,
                "status": match.status,
                "candidate_name": match.candidate_name,
                "assigned_taxid": match.assigned_taxid,
                "assigned_rank": match.assigned_rank,
                "assigned_name": match.assigned_name,
                "hierarchy_match": match.hierarchy_match,
                "hierarchy_family": match.hierarchy_family,
                "hierarchy_order": match.hierarchy_order,
                "hierarchy_class": match.hierarchy_class,
                "hierarchy_superclass": match.hierarchy_superclass,
                "custom_label": match.custom_label,
                "extraction_method": match.extraction_method,
                "original_header": original_header,
                "new_header": new_header,
            }

            write_report_row(report_writer, values)
            write_reference_map_row(refmap_writer, values)

    print(f"Input records: {total}")
    print(f"Already had taxid: {already_had_taxid}")
    print(f"Mapped exact name: {mapped_exact}")
    print(f"Mapped by binomial fallback: {mapped_binomial}")
    print(f"Mapped by genus fallback: {mapped_genus}")
    print(f"Mapped by taxonomy hierarchy fallback: {mapped_hierarchy}")
    print(f"Unresolved/no usable ancestor: {unresolved}")
    print(f"Wrote FASTA: {output_fasta}")
    print(f"Wrote report: {report_tsv}")
    print(f"Wrote reference map: {reference_map_tsv}")


def parse_fallback_ranks(value: str) -> list[str]:
    """Parse comma-delimited fallback ranks."""
    ranks = [rank.strip() for rank in value.split(",") if rank.strip()]
    return ranks or DEFAULT_HIERARCHY_FALLBACK_RANKS


def default_reference_map_path(output_fasta: Path) -> Path:
    """Default reference mapping path beside the rewritten FASTA."""
    return output_fasta.parent / "reference_taxid_mapping.tsv"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Append terminal |taxid= annotations to custom FASTA headers while preserving "
            "original NBDL sequence IDs. Also writes a downstream reference mapping TSV."
        )
    )
    parser.add_argument("-i", "--input-fasta", required=True, type=Path, help="Input FASTA file.")
    parser.add_argument(
        "-o",
        "--output-fasta",
        required=True,
        type=Path,
        help="Output FASTA file with original headers preserved and |taxid= appended where possible.",
    )
    parser.add_argument("-n", "--names2id", required=True, type=Path, help="Two-column TSV: scientific_name<TAB>taxid.")
    parser.add_argument(
        "-r",
        "--report-tsv",
        required=True,
        type=Path,
        help="Output TSV report of exact, fallback, hierarchy, and unresolved records.",
    )
    parser.add_argument(
        "--reference-map",
        default=None,
        type=Path,
        help=(
            "Output downstream reference mapping TSV. Default: reference_taxid_mapping.tsv "
            "beside the output FASTA."
        ),
    )
    parser.add_argument(
        "--taxonomy-hierarchy",
        default=None,
        type=Path,
        help=(
            "Optional fish taxonomy hierarchy CSV. Expected columns include Species, Genus, "
            "Subfamily, Family, Order, Class, SuperClass. Used only after exact/binomial/genus "
            "taxid lookup fails."
        ),
    )
    parser.add_argument(
        "--hierarchy-fallback-ranks",
        default=",".join(DEFAULT_HIERARCHY_FALLBACK_RANKS),
        help=(
            "Comma-delimited hierarchy ranks to try when exact/binomial/genus lookup fail. "
            "Default: Family,Order,Class,SuperClass"
        ),
    )
    parser.add_argument(
        "--replace-existing",
        action="store_true",
        help="Replace existing terminal |taxid= values instead of preserving them.",
    )
    parser.add_argument(
        "--custom-taxon-prefix",
        default="CUSTOM_",
        help=(
            "Prefix for unresolved custom taxa in the report/reference map. These labels are "
            "not appended as |taxid=. Default: CUSTOM_"
        ),
    )

    args = parser.parse_args()

    names2id = load_names2id(args.names2id)
    taxonomy_hierarchy = load_taxonomy_hierarchy(args.taxonomy_hierarchy)
    hierarchy_fallback_ranks = parse_fallback_ranks(args.hierarchy_fallback_ranks)
    reference_map = args.reference_map or default_reference_map_path(args.output_fasta)

    process_fasta(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        names2id=names2id,
        report_tsv=args.report_tsv,
        reference_map_tsv=reference_map,
        taxonomy_hierarchy=taxonomy_hierarchy,
        hierarchy_fallback_ranks=hierarchy_fallback_ranks,
        replace_existing=args.replace_existing,
        custom_taxon_prefix=args.custom_taxon_prefix,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
