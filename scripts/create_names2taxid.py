#!/usr/bin/env python3
"""
Create a names2taxid lookup table from an NCBI/TaxonKit taxdump.

The output is a two-column TSV:

    taxon_name<TAB>taxid

By default, only NCBI scientific names are written. This is the safest lookup
for curated database header manipulation because synonyms and common names can
map ambiguously to multiple taxids.

Examples
--------
Use the default TaxonKit taxdump location:

    python scripts/create_names2taxid.py --output names2taxid.tsv

Use an explicit taxdump directory:

    python scripts/create_names2taxid.py \
        --taxdump ~/.taxonkit \
        --output names2taxid.tsv

Include non-scientific names as well, useful for diagnostics but not generally
recommended for header rewriting:

    python scripts/create_names2taxid.py \
        --include-synonyms \
        --output names2taxid.with_synonyms.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, TextIO


DEFAULT_TAXDUMP_DIR = Path.home() / ".taxonkit"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a two-column names2taxid TSV from an NCBI names.dmp file "
            "as installed/used by TaxonKit."
        )
    )
    parser.add_argument(
        "--taxdump",
        type=Path,
        default=DEFAULT_TAXDUMP_DIR,
        help=(
            "Directory containing names.dmp. Defaults to ~/.taxonkit, which is "
            "the usual TaxonKit taxonomy directory."
        ),
    )
    parser.add_argument(
        "--names-dmp",
        type=Path,
        default=None,
        help=(
            "Optional explicit path to names.dmp. If provided, this overrides "
            "--taxdump."
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("names2taxid.tsv"),
        help="Output TSV path. Defaults to ./names2taxid.tsv.",
    )
    parser.add_argument(
        "--include-synonyms",
        action="store_true",
        help=(
            "Include all name classes, not only scientific names. This can "
            "create ambiguous name-to-taxid mappings and is not recommended "
            "for curated FASTA header rewriting unless you know you need it."
        ),
    )
    parser.add_argument(
        "--duplicate-report",
        type=Path,
        default=None,
        help=(
            "Optional path to write names that map to more than one taxid. "
            "Useful when --include-synonyms is used or when auditing a taxdump."
        ),
    )
    parser.add_argument(
        "--fail-on-duplicates",
        action="store_true",
        help=(
            "Exit with an error if any taxon name maps to more than one taxid. "
            "This is useful for strict database-build workflows."
        ),
    )
    return parser.parse_args()


def split_names_dmp_line(line: str) -> list[str]:
    """Split one NCBI names.dmp line into stripped fields."""
    return [field.strip() for field in line.rstrip("\n").split("|")]


def iter_name_records(names_handle: TextIO) -> Iterable[tuple[str, str, str]]:
    """
    Yield (taxid, name, name_class) records from names.dmp.

    NCBI names.dmp rows have pipe-delimited fields similar to:

        tax_id | name_txt | unique name | name class |
    """
    for line_number, line in enumerate(names_handle, start=1):
        if not line.strip():
            continue

        fields = split_names_dmp_line(line)
        if len(fields) < 4:
            raise ValueError(
                f"Malformed names.dmp line {line_number}: expected at least 4 "
                f"pipe-delimited fields, observed {len(fields)}"
            )

        taxid, name_txt, _unique_name, name_class = fields[:4]
        if not taxid or not name_txt:
            continue

        yield taxid, name_txt, name_class


def build_names2taxid(
    names_dmp: Path,
    include_synonyms: bool = False,
) -> tuple[dict[str, str], dict[str, set[str]], Counter[str]]:
    """
    Build a taxon-name to taxid mapping from names.dmp.

    Returns
    -------
    name_to_taxid:
        Names with a single selected taxid.
    duplicates:
        Names that map to more than one taxid.
    name_class_counts:
        Counts of observed NCBI name classes.
    """
    candidate_taxids: dict[str, set[str]] = defaultdict(set)
    name_class_counts: Counter[str] = Counter()

    with names_dmp.open("r", encoding="utf-8") as handle:
        for taxid, name, name_class in iter_name_records(handle):
            name_class_counts[name_class] += 1

            if not include_synonyms and name_class != "scientific name":
                continue

            candidate_taxids[name].add(taxid)

    duplicates = {
        name: taxids
        for name, taxids in candidate_taxids.items()
        if len(taxids) > 1
    }

    name_to_taxid = {
        name: sorted(taxids, key=lambda value: int(value) if value.isdigit() else value)[0]
        for name, taxids in candidate_taxids.items()
    }

    return name_to_taxid, duplicates, name_class_counts


def write_names2taxid(output: Path, name_to_taxid: dict[str, str]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for name in sorted(name_to_taxid):
            writer.writerow([name, name_to_taxid[name]])


def write_duplicate_report(output: Path, duplicates: dict[str, set[str]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["name", "taxids"])
        for name in sorted(duplicates):
            writer.writerow([name, ";".join(sorted(duplicates[name]))])


def main() -> int:
    args = parse_args()
    names_dmp = args.names_dmp if args.names_dmp else args.taxdump / "names.dmp"

    if not names_dmp.exists():
        print(
            f"ERROR: names.dmp was not found at {names_dmp}\n"
            "Provide --taxdump /path/to/taxdump or --names-dmp /path/to/names.dmp.",
            file=sys.stderr,
        )
        return 1

    name_to_taxid, duplicates, name_class_counts = build_names2taxid(
        names_dmp=names_dmp,
        include_synonyms=args.include_synonyms,
    )

    if duplicates and args.duplicate_report:
        write_duplicate_report(args.duplicate_report, duplicates)

    if duplicates and args.fail_on_duplicates:
        print(
            f"ERROR: {len(duplicates)} names map to more than one taxid. "
            "Rerun with --duplicate-report duplicates.tsv to inspect them.",
            file=sys.stderr,
        )
        return 2

    write_names2taxid(args.output, name_to_taxid)

    selected_mode = "all NCBI name classes" if args.include_synonyms else "scientific names only"
    print(f"Read: {names_dmp}")
    print(f"Wrote: {args.output}")
    print(f"Mode: {selected_mode}")
    print(f"Names written: {len(name_to_taxid)}")

    if duplicates:
        print(f"Ambiguous names: {len(duplicates)}")
        if args.duplicate_report:
            print(f"Duplicate report: {args.duplicate_report}")
        else:
            print("Tip: use --duplicate-report duplicates.tsv to inspect ambiguous names.")

    if name_class_counts:
        print("Observed name classes:")
        for name_class, count in name_class_counts.most_common():
            print(f"  {name_class}: {count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
