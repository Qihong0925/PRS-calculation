#!/usr/bin/env python3
"""Build a PRS samples table by scanning alignment directories.

Directories are given in PRIORITY order: when the same sample appears in more than one
directory, the copy from the earliest directory is used and the others are logged.
Every file found is reported (used / duplicate_ignored); nothing is dropped silently.

Outputs:
  <out>                 sample <TAB> cram <TAB> source_dir   (column 'cram' accepts BAM or CRAM)
  <out>.inventory.tsv   every alignment file found, with status and index check

Usage:
  make_samples_table.py --dirs DIR1 [DIR2 ...] --out samples.tsv
                        [--suffix .dedup.bqsr.bam] [--prefix Homo_sapiens.] [--expect N]
"""
import argparse
import os
import sys

import pandas as pd


def index_for(path):
    for cand in (path + ".bai", path + ".crai", path[:-4] + ".bai", path[:-5] + ".crai"):
        if os.path.exists(cand):
            return cand
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", required=True, help="directories in priority order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--suffix", default=".dedup.bqsr.bam")
    ap.add_argument("--prefix", default="Homo_sapiens.")
    ap.add_argument("--expect", type=int, default=None, help="expected number of unique samples")
    args = ap.parse_args()

    rows = []
    for rank, d in enumerate(args.dirs):
        if not os.path.isdir(d):
            sys.exit(f"Not a directory: {d}")
        files = sorted(f for f in os.listdir(d) if f.endswith(args.suffix))
        if not files:
            print(f"WARNING: no *{args.suffix} in {d}", file=sys.stderr)
        for f in files:
            sample = f[: -len(args.suffix)]
            if args.prefix and sample.startswith(args.prefix):
                sample = sample[len(args.prefix):]
            p = os.path.join(os.path.abspath(d), f)
            rows.append(dict(sample=sample, cram=p, source_dir=os.path.basename(os.path.normpath(d)),
                             priority=rank, size_GB=round(os.path.getsize(p) / 1e9, 2), index=index_for(p)))

    inv = pd.DataFrame(rows).sort_values(["sample", "priority"]).reset_index(drop=True)
    inv["status"] = "used"
    inv.loc[inv.duplicated("sample", keep="first"), "status"] = "duplicate_ignored"
    used = inv[inv.status == "used"]
    dup_samples = inv.loc[inv.status == "duplicate_ignored", "sample"].unique()

    inv.drop(columns="priority").to_csv(args.out + ".inventory.tsv", sep="\t", index=False)
    used[["sample", "cram", "source_dir"]].to_csv(args.out, sep="\t", index=False)

    print("Files found per directory:")
    print(inv.groupby("source_dir", sort=False).size().to_string())
    print(f"\nTotal files: {len(inv)}; unique samples: {len(used)}; duplicated samples: {len(dup_samples)}")
    if len(dup_samples):
        print("\nDuplicates (first listed = used):")
        print(inv[inv["sample"].isin(dup_samples)][["sample", "source_dir", "size_GB", "status"]]
              .to_string(index=False))
    no_idx = used[used["index"] == ""]
    if len(no_idx):
        print(f"\nWARNING: {len(no_idx)} used files have no index:")
        print(no_idx[["sample", "cram"]].to_string(index=False))
    print(f"\nWritten: {args.out}  and  {args.out}.inventory.tsv")
    if args.expect is not None and len(used) != args.expect:
        sys.exit(f"ERROR: expected {args.expect} unique samples, found {len(used)}")


if __name__ == "__main__":
    main()
