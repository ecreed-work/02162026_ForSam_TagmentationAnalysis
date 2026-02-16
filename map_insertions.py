#!/usr/bin/env python3
"""
CHANGE FILE PATHS, COPY AND PASTE IN TERMINAL

Single-sample mode:
python map_insertions.py \
    --fastq <your_fastq> \
    --ref <your_ref> \
    --outdir <your_outdir> \
    --min-mapq 30 \
    --donor-seq <your_donor_seq>

Bulk mode (Excel/CSV):
python map_insertions.py \
    --xlsx /Users/ecreed/Desktop/KelloggRotation/forSam_Tagmentation_Analysis/map_insertions_fastqPath.csv \
    --sheet Sheet1 \
    --min-mapq 30

Bulk mode (auto-detect FASTQ in folder):
python map_insertions.py \
    --fastq-dir /path/to/fastqs \
    --fastq-pattern "*_S*_L001_merged*.fastq*" \
    --ref <your_ref> \
    --outdir <your_outdir> \
    --donor-seq <your_donor_seq> \
    --donor-side 5p \
    --min-mapq 30
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path
from typing import Optional, Tuple, Set

import pandas as pd
import pysam


def run(cmd, *, shell=False):
    cmd_str = cmd if isinstance(cmd, str) else " ".join(map(str, cmd))
    print(f"[cmd] {cmd_str}")
    subprocess.run(cmd, check=True, shell=shell)


def ensure_bwa_index(ref_fa: Path):
    needed = [ref_fa.with_suffix(ref_fa.suffix + ext) for ext in [".amb", ".ann", ".bwt", ".pac", ".sa"]]
    if all(p.exists() for p in needed):
        return
    run(["bwa", "index", str(ref_fa)])


def cutadapt_trim(
    fastq_in: list[Path],
    fastq_out: Path,
    donor_seq: str,
    donor_side: str,
    min_len: int,
    error_rate: float,
    cores: int,
):
    if donor_side == "5p":
        adapter_arg = ["-g", donor_seq]
    else:
        adapter_arg = ["-a", donor_seq]

    cmd = [
        "cutadapt",
        *adapter_arg,
        "--discard-untrimmed",
        "-e", str(error_rate),
        "-m", str(min_len),
        "-j", str(cores),
        "-o", str(fastq_out),
        *map(str, fastq_in),
    ]
    run(cmd)


def bwa_map(ref_fa: Path, fastq: Path, bam_out: Path, cores: int):
    sam_tmp = bam_out.with_suffix(".sam")
    run(["bwa", "mem", "-t", str(cores), str(ref_fa), str(fastq)], shell=False)
    run(f"bwa mem -t {cores} {ref_fa} {fastq} > {sam_tmp}", shell=True)

    run(["samtools", "sort", "-@", str(cores), "-o", str(bam_out), str(sam_tmp)])
    run(["samtools", "index", str(bam_out)])
    try:
        sam_tmp.unlink()
    except FileNotFoundError:
        pass


def ref_aln_span(read: pysam.AlignedSegment) -> Tuple[int, int]:
    return read.reference_start, read.reference_end


def insertion_coord(read: pysam.AlignedSegment, donor_side: str) -> Optional[int]:
    if read.is_unmapped:
        return None
    if read.reference_start is None or read.reference_end is None:
        return None
    if read.cigarstring is None or read.cigarstring == "*":
        return None

    start0, end0_excl = ref_aln_span(read)
    if end0_excl <= start0:
        return None

    forward = not read.is_reverse
    return start0 if (donor_side == "5p" and forward) else (
        end0_excl - 1 if donor_side == "5p" else (end0_excl - 1 if forward else start0)
    )


def call_insertions(
    bam_path: Path,
    out_tsv: Path,
    donor_side: str,
    min_mapq: int,
    primary_only: bool,
):
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    with open(out_tsv, "w") as f:
        f.write("\t".join([
            "read_id", "ref", "ins0", "ins1", "strand", "mapq", "cigar",
            "aln_start0", "aln_end0_excl"
        ]) + "\n")

        n_total = 0
        n_kept = 0

        for read in bam.fetch(until_eof=True):
            n_total += 1
            if read.is_unmapped:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if primary_only and (read.is_secondary or read.is_supplementary):
                continue

            ins0 = insertion_coord(read, donor_side=donor_side)
            if ins0 is None:
                continue

            ref = bam.get_reference_name(read.reference_id)
            strand = "-" if read.is_reverse else "+"
            start0, end0_excl = ref_aln_span(read)

            f.write("\t".join(map(str, [
                read.query_name, ref, ins0, ins0 + 1, strand, read.mapping_quality,
                read.cigarstring, start0, end0_excl
            ])) + "\n")
            n_kept += 1

    bam.close()
    print(f"[info] BAM reads scanned: {n_total}")
    print(f"[info] Insertions written: {n_kept}")
    print(f"[info] Output TSV: {out_tsv}")


def run_one_sample(
    *,
    fastq: list[Path],
    ref: Path,
    outdir: Path,
    donor_seq: str,
    donor_side: str,
    min_len: int,
    error_rate: float,
    cores: int,
    min_mapq: int,
    primary_only: bool,
    indexed_refs: Set[Path],
):
    outdir.mkdir(parents=True, exist_ok=True)

    fastq = [p.expanduser().resolve() for p in fastq]
    ref = ref.expanduser().resolve()

    stem = fastq[0].stem
    trimmed_fastq = outdir / f"{stem}.trimmed.fastq"
    bam = outdir / f"{stem}.mapped.sorted.bam"
    out_tsv = outdir / f"{stem}.insertions.{donor_side}.tsv"

    if ref not in indexed_refs:
        ensure_bwa_index(ref)
        indexed_refs.add(ref)

    cutadapt_trim(
        fastq_in=fastq,
        fastq_out=trimmed_fastq,
        donor_seq=donor_seq,
        donor_side=donor_side,
        min_len=min_len,
        error_rate=error_rate,
        cores=cores,
    )

    bwa_map(ref_fa=ref, fastq=trimmed_fastq, bam_out=bam, cores=cores)

    call_insertions(
        bam_path=bam,
        out_tsv=out_tsv,
        donor_side=donor_side,
        min_mapq=min_mapq,
        primary_only=primary_only,
    )


def _parse_fastq_list(value: str) -> list[Path]:
    parts = [p.strip() for p in re.split(r"[;,]+", value.strip()) if p.strip()]
    return [Path(p) for p in parts]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", type=Path)
    ap.add_argument("--sheet", default="Sheet1")
    ap.add_argument("--fastq", type=Path)
    ap.add_argument("--fastq-dir", type=Path)
    ap.add_argument("--fastq-pattern", default="*_S*_L001_merged*.fastq*")
    ap.add_argument("--ref", type=Path)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--donor-seq", dest="donor_seq")
    ap.add_argument("--donor-side", dest="donor_side", choices=["5p", "3p"], default="5p")
    ap.add_argument("--min-len", type=int, default=20)
    ap.add_argument("--error-rate", type=float, default=0.1)
    ap.add_argument("--cores", type=int, default=8)
    ap.add_argument("--min-mapq", type=int, default=20)
    ap.add_argument("--primary-only", action="store_true")

    args = ap.parse_args()
    indexed_refs: Set[Path] = set()

    if args.xlsx is not None:
        if args.xlsx.suffix.lower() == ".csv":
            df = pd.read_csv(args.xlsx)
        else:
            df = pd.read_excel(args.xlsx, sheet_name=args.sheet)

        for i, row in df.iterrows():
            if row.isna().all():
                continue
            fastq_raw = str(row["fastq"]).strip() if not pd.isna(row["fastq"]) else None
            ref = Path(str(row["ref"]).strip()) if not pd.isna(row["ref"]) else None
            outdir = Path(str(row["outdir"]).strip()) if not pd.isna(row["outdir"]) else None
            donor_seq = str(row["donor-seq"]).strip() if not pd.isna(row["donor-seq"]) else None
            donor_side = str(row["donor-side"]).strip() if not pd.isna(row["donor-side"]) else "5p"

            if fastq_raw is None or ref is None or outdir is None or donor_seq is None:
                raise SystemExit(f"[error] Row {i+2} has missing required values.")

            fastq_list = _parse_fastq_list(fastq_raw)

            run_one_sample(
                fastq=fastq_list,
                ref=ref,
                outdir=outdir,
                donor_seq=donor_seq,
                donor_side=donor_side,
                min_len=args.min_len,
                error_rate=args.error_rate,
                cores=args.cores,
                min_mapq=args.min_mapq,
                primary_only=args.primary_only,
                indexed_refs=indexed_refs,
            )
        return

    raise SystemExit("[error] Use --xlsx for bulk mode.")


if __name__ == "__main__":
    main()
