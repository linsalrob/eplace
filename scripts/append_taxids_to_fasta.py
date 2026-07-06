#!/usr/bin/env python3

"""
Append usable NCBI taxids to custom reference FASTA headers and write a BLAST taxid map.

This utility is intended for custom reference databases such as CSIRO/NBDL,
BOLD, institutional FASTA exports, or mixed curated databases where sequence
headers contain organism names but may not already contain NCBI taxids.

Resolution ladder:

1. Extract organism/species name from flexible FASTA header formats.
2. Map the exact name to a taxid using names2id.tsv.
3. If exact name is absent, fall back to genus-level taxid when genus exists.
4. If genus is absent and --taxonomy-hierarchy is supplied, use that hierarchy
   to find family/order/class/superclass for the species or genus, then map the
   first available ancestor name to names2id.tsv.
5. If no valid NCBI ancestor can be found, preserve the custom organism name and
   assign a stable custom label in metadata/report, but do not append fake
   |taxid= values.

The script writes three outputs:

1. Rewritten FASTA with ePLACE metadata and usable |taxid= annotations where possible.
2. Report TSV describing exact, fallback, and unresolved mappings.
3. BLAST taxid map TSV, defaulting to taxid_map.tsv beside the output FASTA.

Expected names2id.tsv format:

    scientific_name<TAB>taxid

Expected --taxonomy-hierarchy CSV columns, case-insensitive:

    Species, Genus, Subfamily, Family, Order, Class, SuperClass

Extra columns are allowed. The script uses Species and Genus as lookup keys and
tries fallback ranks in this order by default:

    Family -> Order -> Class -> SuperClass

Important:

- The appended |taxid= value should always be a real taxid present in your
  taxonomy dump. Do not append fake custom IDs as |taxid= unless your taxonomy
  database has been extended to include them.
- For custom species absent from NCBI, the script preserves the original name in
  [eplace_original_organism=...] even when |taxid= points to a conservative
  genus/family/order/class/superclass ancestor.
- The taxid_map.tsv file is what makeblastdb needs for -taxid_map. The |taxid=
  annotation in the FASTA header is useful metadata, but BLAST does not reliably
  convert arbitrary header text into internal sequence taxids without a taxid map.
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

HEADER_METADATA_KEYS = [
    "eplace_original_organism",
    "eplace_taxid_assignment",
    "eplace_taxid_rank",
    "eplace_taxid_name",
    "eplace_hierarchy_match",
    "eplace_custom_label",
]

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
    custom_label: str = ""


def normalize_column_name(value: str) -> str:
    """Normalize CSV column names for case-insensitive matching."""
    return re.sub(r"[^a-z0-9]", "", value.lower())


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


def row_value(row: dict[str, str], column_map: dict[str, str], desired_column: str) -> str:
    """Return a CSV row value by case-insensitive normalized column name."""
    actual = column_map.get(normalize_column_name(desired_column))
    if not actual:
        return ""
    return row.get(actual, "").strip()


def load_taxonomy_hierarchy(path: Optional[Path]) -> TaxonomyHierarchy:
    """
    Load an optional fish taxonomy hierarchy CSV.

    Expected useful columns are Species, Genus, Subfamily, Family, Order, Class,
    and SuperClass. Column names are matched case-insensitively.
    """
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

            # Store the first genus-level record only. If a genus appears multiple
            # times, family/order/class should normally be stable; preserving the
            # first avoids later species overwriting earlier genus metadata.
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


def strip_eplace_metadata(header: str) -> str:
    """Remove existing ePLACE taxid-assignment metadata fields."""
    cleaned = header.rstrip()
    for key in HEADER_METADATA_KEYS:
        cleaned = re.sub(rf"\s*\[{re.escape(key)}=[^\]]*\]", "", cleaned)
    return cleaned


def bracket_escape(value: str) -> str:
    """Make a value safe for bracket-style FASTA metadata."""
    return value.replace("]", ")").replace("[", "(").strip()


def sequence_id_from_header(header: str) -> str:
    """
    Return the sequence ID BLAST will see when -parse_seqids is used.

    This is the first whitespace-delimited token after removing the leading >.
    The taxid_map.tsv first column must match this token exactly.
    """
    return header.lstrip(">").strip().split()[0]


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

    if match.hierarchy_match:
        fields.append(f"[eplace_hierarchy_match={bracket_escape(match.hierarchy_match)}]")

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


def try_hierarchy_fallback(
    candidate_name: str,
    names2id: dict[str, str],
    hierarchy: TaxonomyHierarchy,
    fallback_ranks: list[str],
) -> Optional[NameMatch]:
    """
    Try mapping candidate species/genus to an ancestor using taxonomy hierarchy.

    The first ancestor rank whose name exists in names2id.tsv is used.
    """
    record = hierarchy.lookup(candidate_name)
    if record is None:
        return None

    hierarchy_identity = record.species or record.genus or candidate_name

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
            hierarchy_match=hierarchy_identity,
        )

    return None


def resolve_taxid(
    candidate_name: Optional[str],
    names2id: dict[str, str],
    *,
    taxonomy_hierarchy: Optional[TaxonomyHierarchy] = None,
    hierarchy_fallback_ranks: Optional[list[str]] = None,
    custom_taxon_prefix: Optional[str] = None,
    custom_taxa_seen: Optional[dict[str, str]] = None,
) -> NameMatch:
    """
    Resolve candidate name to exact, genus, hierarchy-derived ancestor, or unmapped.

    The fallback behaviour is intentionally conservative: the original species
    name is preserved in ePLACE metadata while |taxid= points to a real NCBI
    ancestor such as genus, family, order, class, or superclass.
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

    if taxonomy_hierarchy is not None:
        hierarchy_match = try_hierarchy_fallback(
            candidate_name=candidate_name,
            names2id=names2id,
            hierarchy=taxonomy_hierarchy,
            fallback_ranks=hierarchy_fallback_ranks or DEFAULT_HIERARCHY_FALLBACK_RANKS,
        )
        if hierarchy_match:
            return hierarchy_match

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
    taxid_map: Path,
    *,
    taxonomy_hierarchy: Optional[TaxonomyHierarchy] = None,
    hierarchy_fallback_ranks: Optional[list[str]] = None,
    replace_existing: bool = False,
    refresh_eplace_metadata: bool = False,
    custom_taxon_prefix: Optional[str] = None,
) -> None:
    """Process FASTA headers and append usable taxids plus ePLACE metadata."""
    total = 0
    already_had_taxid = 0
    mapped_exact = 0
    mapped_genus = 0
    mapped_hierarchy = 0
    unresolved = 0
    taxid_map_written = 0
    custom_taxa_seen: dict[str, str] = {}

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    report_tsv.parent.mkdir(parents=True, exist_ok=True)
    taxid_map.parent.mkdir(parents=True, exist_ok=True)

    with (
        input_fasta.open() as inp,
        output_fasta.open("w") as out,
        report_tsv.open("w") as rep,
        taxid_map.open("w") as taxmap,
    ):
        report_writer = csv.writer(rep, delimiter="\t")
        report_writer.writerow([
            "seq_id",
            "status",
            "candidate_name",
            "assigned_taxid",
            "assigned_rank",
            "assigned_name",
            "hierarchy_match",
            "custom_label",
            "extraction_method",
            "taxid_map_written",
            "original_header",
            "new_header",
        ])

        for line in inp:
            if not line.startswith(">"):
                out.write(line)
                continue

            total += 1
            original_header = line.rstrip("\n")
            seq_id = sequence_id_from_header(original_header)

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
                taxmap.write(f"{seq_id}\t{current_taxid}\n")
                taxid_map_written += 1

                report_writer.writerow([
                    seq_id,
                    match.status,
                    match.candidate_name,
                    match.assigned_taxid,
                    match.assigned_rank,
                    match.assigned_name,
                    match.hierarchy_match,
                    match.custom_label,
                    method,
                    "Yes",
                    original_header,
                    working_header,
                ])
                continue

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
            elif match.status == "mapped_genus_fallback":
                mapped_genus += 1
            elif "hierarchy_fallback" in match.status:
                mapped_hierarchy += 1
            else:
                unresolved += 1

            working_header = append_assignment_metadata(working_header, match)
            working_header = append_taxid_to_header(
                working_header,
                match.assigned_taxid or None,
                replace_existing=replace_existing,
            )

            taxid_map_status = "No"
            if match.assigned_taxid:
                taxmap.write(f"{seq_id}\t{match.assigned_taxid}\n")
                taxid_map_written += 1
                taxid_map_status = "Yes"

            out.write(working_header + "\n")
            report_writer.writerow([
                seq_id,
                match.status,
                match.candidate_name,
                match.assigned_taxid,
                match.assigned_rank,
                match.assigned_name,
                match.hierarchy_match,
                match.custom_label,
                method,
                taxid_map_status,
                original_header,
                working_header,
            ])

    print(f"Input records: {total}")
    print(f"Already had taxid: {already_had_taxid}")
    print(f"Mapped exact name: {mapped_exact}")
    print(f"Mapped by genus fallback: {mapped_genus}")
    print(f"Mapped by taxonomy hierarchy fallback: {mapped_hierarchy}")
    print(f"Unresolved/no usable ancestor: {unresolved}")
    print(f"Taxid map rows written: {taxid_map_written}")
    print(f"Wrote FASTA: {output_fasta}")
    print(f"Wrote report: {report_tsv}")
    print(f"Wrote BLAST taxid map: {taxid_map}")


def parse_fallback_ranks(value: str) -> list[str]:
    """Parse comma-delimited fallback ranks."""
    ranks = [rank.strip() for rank in value.split(",") if rank.strip()]
    return ranks or DEFAULT_HIERARCHY_FALLBACK_RANKS


def default_taxid_map_path(output_fasta: Path) -> Path:
    """Default to taxid_map.tsv beside the rewritten FASTA."""
    return output_fasta.parent / "taxid_map.tsv"


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Append NCBI taxids to custom FASTA headers using organism/species "
            "names, with genus and optional taxonomy-hierarchy fallback. Also "
            "writes a BLAST-compatible taxid_map.tsv."
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
        help="Output TSV report of exact, genus-fallback, hierarchy-fallback, and unmapped records.",
    )
    parser.add_argument(
        "--taxid-map",
        default=None,
        type=Path,
        help=(
            "Output BLAST taxid map TSV. Default: taxid_map.tsv beside the output FASTA. "
            "Use this file with makeblastdb -taxid_map."
        ),
    )
    parser.add_argument(
        "--taxonomy-hierarchy",
        default=None,
        type=Path,
        help=(
            "Optional fish taxonomy hierarchy CSV. Expected columns include "
            "Species, Genus, Subfamily, Family, Order, Class, SuperClass. "
            "Used only after exact species/name and genus taxid lookup fail."
        ),
    )
    parser.add_argument(
        "--hierarchy-fallback-ranks",
        default=",".join(DEFAULT_HIERARCHY_FALLBACK_RANKS),
        help=(
            "Comma-delimited hierarchy ranks to try when exact/genus lookup fail. "
            "Default: Family,Order,Class,SuperClass"
        ),
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
    taxonomy_hierarchy = load_taxonomy_hierarchy(args.taxonomy_hierarchy)
    hierarchy_fallback_ranks = parse_fallback_ranks(args.hierarchy_fallback_ranks)
    taxid_map = args.taxid_map or default_taxid_map_path(args.output_fasta)

    process_fasta(
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta,
        names2id=names2id,
        report_tsv=args.report_tsv,
        taxid_map=taxid_map,
        taxonomy_hierarchy=taxonomy_hierarchy,
        hierarchy_fallback_ranks=hierarchy_fallback_ranks,
        replace_existing=args.replace_existing,
        refresh_eplace_metadata=args.refresh_eplace_metadata,
        custom_taxon_prefix=args.custom_taxon_prefix,
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
