"""Standalone amplicon preprocessing for ePLACE.

Demultiplex paired FASTQ reads using separate I1/I2 index FASTQ files, trim
primers with Cutadapt, infer ASVs with DADA2, and export an ePLACE-ready FASTA
plus an ASV-by-sample count table.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import logging
import shutil
import subprocess
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, TextIO, Tuple

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class SampleIndex:
    sample_id: str
    index1: str
    index2: str


def _open_text(path: Path, mode: str = "rt") -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def _fastq_records(handle: TextIO, source: Path) -> Iterator[Tuple[str, str, str, str]]:
    while True:
        header = handle.readline()
        if not header:
            return
        sequence = handle.readline()
        plus = handle.readline()
        quality = handle.readline()
        if not sequence or not plus or not quality:
            raise ValueError(f"Truncated FASTQ record in {source}")
        if not header.startswith("@") or not plus.startswith("+"):
            raise ValueError(f"Invalid FASTQ record in {source}: {header.rstrip()}")
        if len(sequence.rstrip()) != len(quality.rstrip()):
            raise ValueError(f"Sequence/quality length mismatch in {source}: {header.rstrip()}")
        yield header, sequence, plus, quality


def _read_id(header: str) -> str:
    token = header[1:].split()[0]
    if token.endswith("/1") or token.endswith("/2"):
        token = token[:-2]
    return token


def _hamming(observed: str, expected: str) -> Optional[int]:
    observed = observed.upper()
    expected = expected.upper()
    if len(observed) < len(expected):
        return None
    observed = observed[: len(expected)]
    return sum(a != b for a, b in zip(observed, expected))


def _load_sample_sheet(path: Path) -> List[SampleIndex]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "index1", "index2"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"Index sheet is missing columns: {', '.join(sorted(missing))}. "
                "Required columns are sample_id, index1, index2."
            )
        samples = [
            SampleIndex(
                row["sample_id"].strip(),
                row["index1"].strip().upper(),
                row["index2"].strip().upper(),
            )
            for row in reader
            if row.get("sample_id", "").strip()
        ]
    if not samples:
        raise ValueError("Index sheet contains no samples")
    ids = [sample.sample_id for sample in samples]
    if len(ids) != len(set(ids)):
        raise ValueError("Index sheet contains duplicate sample_id values")
    pairs = [(sample.index1, sample.index2) for sample in samples]
    if len(pairs) != len(set(pairs)):
        raise ValueError("Index sheet contains duplicate index1/index2 combinations")
    return samples


def _match_sample(index1: str, index2: str, samples: List[SampleIndex], mismatches: int) -> Optional[str]:
    candidates: List[Tuple[int, str]] = []
    for sample in samples:
        d1 = _hamming(index1, sample.index1)
        d2 = _hamming(index2, sample.index2)
        if d1 is None or d2 is None or d1 > mismatches or d2 > mismatches:
            continue
        candidates.append((d1 + d2, sample.sample_id))
    if not candidates:
        return None
    candidates.sort()
    if len(candidates) > 1 and candidates[0][0] == candidates[1][0]:
        return None
    return candidates[0][1]


def demultiplex(args: argparse.Namespace, samples: List[SampleIndex], output_dir: Path) -> Path:
    demux_dir = output_dir / "01_demultiplexed"
    demux_dir.mkdir(parents=True, exist_ok=True)
    counts = {sample.sample_id: 0 for sample in samples}
    unmatched = 0
    invalid_read_ids = 0

    with ExitStack() as stack:
        input_paths = (args.reads_r1, args.reads_r2, args.index_i1, args.index_i2)
        inputs = [stack.enter_context(_open_text(path)) for path in input_paths]
        iterators = [_fastq_records(handle, path) for handle, path in zip(inputs, input_paths)]
        writers = {
            sample.sample_id: (
                stack.enter_context(gzip.open(demux_dir / f"{sample.sample_id}_R1.fastq.gz", "wt")),
                stack.enter_context(gzip.open(demux_dir / f"{sample.sample_id}_R2.fastq.gz", "wt")),
            )
            for sample in samples
        }

        while True:
            records = []
            ended = []
            for iterator in iterators:
                try:
                    records.append(next(iterator))
                    ended.append(False)
                except StopIteration:
                    records.append(None)
                    ended.append(True)
            if all(ended):
                break
            if any(ended):
                raise ValueError("R1, R2, I1 and I2 FASTQ files contain different numbers of records")
            r1, r2, i1, i2 = records
            assert r1 and r2 and i1 and i2
            ids = {_read_id(record[0]) for record in (r1, r2, i1, i2)}
            if len(ids) != 1:
                invalid_read_ids += 1
                continue
            sample_id = _match_sample(i1[1].strip(), i2[1].strip(), samples, args.index_mismatches)
            if sample_id is None:
                unmatched += 1
                continue
            writers[sample_id][0].writelines(r1)
            writers[sample_id][1].writelines(r2)
            counts[sample_id] += 1

    report = output_dir / "demultiplexing_summary.tsv"
    with report.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "assigned_read_pairs"])
        for sample in samples:
            writer.writerow([sample.sample_id, counts[sample.sample_id]])
        writer.writerow(["__unmatched_or_tied__", unmatched])
        writer.writerow(["__invalid_read_ids__", invalid_read_ids])
    logger.info("Demultiplexed reads written to %s", demux_dir)
    return demux_dir


def trim_primers(args: argparse.Namespace, samples: List[SampleIndex], demux_dir: Path, output_dir: Path) -> Path:
    trimmed_dir = output_dir / "02_trimmed"
    log_dir = output_dir / "logs" / "cutadapt"
    trimmed_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    for sample in samples:
        input_r1 = demux_dir / f"{sample.sample_id}_R1.fastq.gz"
        input_r2 = demux_dir / f"{sample.sample_id}_R2.fastq.gz"
        output_r1 = trimmed_dir / f"{sample.sample_id}.trimmed_R1.fastq.gz"
        output_r2 = trimmed_dir / f"{sample.sample_id}.trimmed_R2.fastq.gz"
        command = [
            "cutadapt", "-j", str(args.threads),
            "-g", f"^{args.forward_primer}", "-G", f"^{args.reverse_primer}",
            "-e", str(args.primer_error_rate), "--no-indels", "--discard-untrimmed",
            "--pair-filter=any", "--minimum-length", str(args.minimum_length),
            "-o", str(output_r1), "-p", str(output_r2), str(input_r1), str(input_r2),
        ]
        with (log_dir / f"{sample.sample_id}.log").open("wt", encoding="utf-8") as log_handle:
            subprocess.run(command, check=True, stdout=log_handle, stderr=subprocess.STDOUT)
    logger.info("Primer-trimmed reads written to %s", trimmed_dir)
    return trimmed_dir


def run_dada2(args: argparse.Namespace, trimmed_dir: Path, output_dir: Path) -> None:
    r_script = Path(__file__).with_name("dada2_preprocess.R")
    command = [
        "Rscript", str(r_script), "--input-dir", str(trimmed_dir),
        "--output-dir", str(output_dir / "03_dada2"), "--threads", str(args.threads),
        "--trunc-len-f", str(args.trunc_len_f), "--trunc-len-r", str(args.trunc_len_r),
        "--max-ee-f", str(args.max_ee_f), "--max-ee-r", str(args.max_ee_r),
        "--trunc-q", str(args.trunc_q), "--min-length", str(args.minimum_length),
        "--min-overlap", str(args.min_overlap), "--max-mismatch", str(args.max_mismatch),
        "--pool", args.pool,
    ]
    subprocess.run(command, check=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="eplace-preprocess",
        description=(
            "Demultiplex paired amplicon FASTQ reads from I1/I2 indexes, trim primers "
            "with Cutadapt, infer ASVs with DADA2, and export ASV sequences and counts."
        ),
    )
    parser.add_argument("--reads-r1", type=Path, required=True)
    parser.add_argument("--reads-r2", type=Path, required=True)
    parser.add_argument("--index-i1", type=Path, required=True)
    parser.add_argument("--index-i2", type=Path, required=True)
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--forward-primer", required=True)
    parser.add_argument("--reverse-primer", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--index-mismatches", type=int, default=0)
    parser.add_argument("--primer-error-rate", type=float, default=0.15)
    parser.add_argument("--minimum-length", type=int, default=50)
    parser.add_argument("--trunc-len-f", type=int, default=0)
    parser.add_argument("--trunc-len-r", type=int, default=0)
    parser.add_argument("--max-ee-f", type=float, default=2.0)
    parser.add_argument("--max-ee-r", type=float, default=6.0)
    parser.add_argument("--trunc-q", type=int, default=6)
    parser.add_argument("--min-overlap", type=int, default=10)
    parser.add_argument("--max-mismatch", type=int, default=1)
    parser.add_argument("--pool", choices=("false", "pseudo", "true"), default="pseudo")
    return parser


def main(argv: Optional[List[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)
    for path in (args.reads_r1, args.reads_r2, args.index_i1, args.index_i2, args.sample_sheet):
        if not path.is_file():
            parser.error(f"Input file does not exist: {path}")
    if args.threads < 1:
        parser.error("--threads must be at least 1")
    if args.index_mismatches < 0:
        parser.error("--index-mismatches cannot be negative")
    for executable in ("cutadapt", "Rscript"):
        if shutil.which(executable) is None:
            parser.error(f"Required executable not found on PATH: {executable}")

    try:
        samples = _load_sample_sheet(args.sample_sheet)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        demux_dir = demultiplex(args, samples, args.output_dir)
        trimmed_dir = trim_primers(args, samples, demux_dir, args.output_dir)
        run_dada2(args, trimmed_dir, args.output_dir)
    except (OSError, ValueError, subprocess.CalledProcessError) as exc:
        logger.error("Preprocessing failed: %s", exc)
        return 1

    logger.info("Preprocessing complete. Final outputs are in %s", args.output_dir / "03_dada2")
    return 0


if __name__ == "__main__":
    sys.exit(main())
