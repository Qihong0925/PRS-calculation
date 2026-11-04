#!/usr/bin/env python3
"""Build GRCh38 PRS sites and a plink2 score file from a CanRisk PRS file (GRCh37).

Steps:
  1. Parse CanRisk file: line 1 'alpha = x', then CSV with explicit
     Reference_Allele (non-effect allele) and Effect_Allele.
  2. Lift GRCh37 -> GRCh38 with the UCSC chain (pyliftover). SNVs landing on the
     '-' strand are reverse-complemented. Indels on the '-' strand cannot be lifted
     by coordinate alone and require a documented manual correction.
  3. Apply manual corrections (keyed on GRCh37 chr/pos/ref/effect).
  4. Check REF against the GRCh38 FASTA. If the GRCh38 reference carries the effect
     allele (reference changed between builds), REF/ALT are swapped and the effect
     allele is kept.
  5. Left-align/trim with bcftools norm.
  6. Cross-check GRCh38 position against the PGS Catalog harmonised file (optional).
  7. site_status.tsv covers EVERY input variant; nothing is dropped silently.
"""
import argparse
import os
import subprocess
import sys
import tempfile

import pandas as pd
from pyliftover import LiftOver

CANRISK_COLS = ["Chromosome", "Position", "Reference_Allele", "Effect_Allele",
                "Effect_Allele_Frequency", "Log_Odds_Ratio"]
CORR_COLS = ["chr37", "pos37", "ref37", "effect37", "chrom38", "pos38",
             "ref38", "alt38", "effect38", "rsid", "evidence"]
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_canrisk(path):
    with open(path) as fh:
        first = fh.readline().strip()
        if not first.lower().startswith("alpha"):
            sys.exit(f"Unexpected first line in {path}: {first}")
        alpha = float(first.split("=")[1])
        df = pd.read_csv(fh, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    missing = [c for c in CANRISK_COLS if c not in df.columns]
    if missing:
        sys.exit(f"CanRisk file missing columns: {missing}")
    for c in CANRISK_COLS:
        df[c] = df[c].str.strip()
    df["Reference_Allele"] = df["Reference_Allele"].str.upper()
    df["Effect_Allele"] = df["Effect_Allele"].str.upper()
    return df, alpha


def read_corrections(path):
    if not path or path == "none":
        return {}
    corr = pd.read_csv(path, sep="\t", dtype=str)
    missing = [c for c in CORR_COLS if c not in corr.columns]
    if missing:
        sys.exit(f"Corrections file missing columns: {missing}")
    out = {}
    for _, c in corr.iterrows():
        key = (c["chr37"], c["pos37"], c["ref37"].upper(), c["effect37"].upper())
        out[key] = c.to_dict()
    return out


def read_pgs_positions(path):
    if not path or path == "none" or not os.path.exists(path):
        return None
    pgs = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    return {(r["chr_name"], r["chr_position"]): r["hm_pos"] for _, r in pgs.iterrows()}


def read_fai(fasta):
    contigs = []
    with open(fasta + ".fai") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            contigs.append((f[0], int(f[1])))
    return contigs


def fetch_ref_seqs(fasta, regions, tmpdir):
    """Return reference sequences (upper case) for regions, in input order."""
    region_file = os.path.join(tmpdir, "regions.txt")
    with open(region_file, "w") as fh:
        fh.write("\n".join(regions) + "\n")
    out = subprocess.run(["samtools", "faidx", fasta, "-r", region_file],
                         capture_output=True, text=True, check=True).stdout
    seqs, cur = [], None
    for line in out.splitlines():
        if line.startswith(">"):
            if cur is not None:
                seqs.append(cur)
            cur = ""
        else:
            cur += line.strip().upper()
    if cur is not None:
        seqs.append(cur)
    if len(seqs) != len(regions):
        sys.exit(f"faidx returned {len(seqs)} sequences for {len(regions)} regions")
    return seqs


def bcftools_norm(raw_vcf, fasta, tmpdir):
    """Left-align/trim with bcftools norm; return {ID: (chrom, pos, ref, alt)}."""
    cmd = (f"bcftools norm -f {fasta} -c e -Ou {raw_vcf} "
           f"| bcftools sort -T {tmpdir} -Ou "
           f"| bcftools query -f '%ID\\t%CHROM\\t%POS\\t%REF\\t%ALT\\n'")
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if p.returncode != 0:
        sys.exit(f"bcftools norm failed:\n{p.stderr}")
    print("bcftools norm:", " ".join(l for l in p.stderr.splitlines() if l.startswith("Lines")))
    out = {}
    for line in p.stdout.splitlines():
        vid, chrom, pos, ref, alt = line.split("\t")
        out[vid] = (chrom, int(pos), ref, alt)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--canrisk", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--ref-fasta", required=True)
    ap.add_argument("--corrections", default="none")
    ap.add_argument("--pgs-check", default="none")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    df, alpha = read_canrisk(args.canrisk)
    df["row_id"] = [f"v{i:04d}" for i in range(len(df))]
    corr = read_corrections(args.corrections)
    pgs_pos = read_pgs_positions(args.pgs_check)
    lo = LiftOver(args.chain)

    contigs = read_fai(args.ref_fasta)
    contig_rank = {c: i for i, (c, _) in enumerate(contigs)}
    contig_len = dict(contigs)
    if "chr1" not in contig_rank:
        sys.exit("Reference FASTA has no 'chr1'; this script expects chr-prefixed GRCh38 contigs")

    status, cand, used_corr, pgs_note = {}, {}, set(), {}
    for _, r in df.iterrows():
        rid = r["row_id"]
        key = (r["Chromosome"], r["Position"], r["Reference_Allele"], r["Effect_Allele"])

        if key in corr:
            c = corr[key]
            used_corr.add(key)
            cand[rid] = dict(chrom=c["chrom38"], pos=int(c["pos38"]), ref=c["ref38"].upper(),
                             alt=c["alt38"].upper(), ea=c["effect38"].upper(),
                             method="manual_correction", manual=True)
        else:
            chrom37 = "chr" + r["Chromosome"]
            hits = lo.convert_coordinate(chrom37, int(r["Position"]) - 1)
            if not hits:
                status[rid] = ("unresolved", "liftover_failed"); continue
            if len(hits) > 1:
                status[rid] = ("unresolved", f"liftover_multiple_hits:{len(hits)}"); continue
            chrom, pos0, strand, _ = hits[0]
            if chrom != chrom37:
                status[rid] = ("unresolved", f"liftover_other_contig:{chrom}"); continue
            if chrom not in contig_rank:
                status[rid] = ("unresolved", f"contig_not_in_reference:{chrom}"); continue
            ref, eff = r["Reference_Allele"], r["Effect_Allele"]
            if strand == "-":
                if len(ref) == 1 and len(eff) == 1:
                    ref, eff, method = revcomp(ref), revcomp(eff), "liftover_minus_strand_revcomp"
                else:
                    status[rid] = ("unresolved", "indel_on_minus_strand_needs_manual_correction"); continue
            else:
                method = "liftover"
            cand[rid] = dict(chrom=chrom, pos=pos0 + 1, ref=ref, alt=eff, ea=eff,
                             method=method, manual=False)

        if pgs_pos is not None:
            hm = pgs_pos.get((r["Chromosome"], r["Position"]))
            if hm is None or pd.isna(hm):
                pgs_note[rid] = "not_in_pgs"
            else:
                pgs_note[rid] = "same_pos" if int(hm) == cand[rid]["pos"] else f"pgs_pos_differs:{hm}"

    unused = set(corr) - used_corr
    if unused:
        sys.exit(f"Manual corrections did not match any CanRisk record: {sorted(unused)}")

    resolved = {}
    with tempfile.TemporaryDirectory(dir=args.outdir) as tmpdir:
        items = list(cand.items())
        regions = [f"{v['chrom']}:{v['pos']}-{v['pos'] + max(len(v['ref']), len(v['alt'])) - 1}"
                   for _, v in items]
        seqs = fetch_ref_seqs(args.ref_fasta, regions, tmpdir) if regions else []

        for (rid, v), seq in zip(items, seqs):
            if seq.startswith(v["ref"]):
                ref, alt, note = v["ref"], v["alt"], "ref_match"
            elif not v["manual"] and seq.startswith(v["alt"]):
                ref, alt, note = v["alt"], v["ref"], "grch38_ref_is_effect_allele"
            else:
                status[rid] = ("unresolved", f"ref_mismatch:ref_seq={seq}")
                continue
            resolved[rid] = dict(chrom=v["chrom"], pos=v["pos"], ref=ref, alt=alt,
                                 ea_is_ref=(v["ea"] == ref), notes=[v["method"], note])

        raw_vcf = os.path.join(tmpdir, "raw.vcf")
        used = sorted({v["chrom"] for v in resolved.values()}, key=contig_rank.get)
        with open(raw_vcf, "w") as fh:
            fh.write("##fileformat=VCFv4.2\n")
            for c in used:
                fh.write(f"##contig=<ID={c},length={contig_len[c]}>\n")
            fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for rid, v in sorted(resolved.items(), key=lambda kv: (contig_rank[kv[1]["chrom"]], kv[1]["pos"])):
                fh.write(f"{v['chrom']}\t{v['pos']}\t{rid}\t{v['ref']}\t{v['alt']}\t.\t.\t.\n")
        normed = bcftools_norm(raw_vcf, args.ref_fasta, tmpdir) if resolved else {}

    lost = set(resolved) - set(normed)
    if lost:
        sys.exit(f"bcftools norm dropped {len(lost)} variants: {sorted(lost)[:10]}")

    info = df.set_index("row_id")
    sites, scores = [], []
    for rid, v in resolved.items():
        chrom, pos, ref, alt = normed[rid]
        ea = ref if v["ea_is_ref"] else alt
        oa = alt if v["ea_is_ref"] else ref
        notes = list(v["notes"])
        if (pos, ref, alt) != (v["pos"], v["ref"], v["alt"]):
            notes.append(f"normalised_from:{v['pos']}:{v['ref']}:{v['alt']}")
        vid = f"{chrom}:{pos}:{ref}:{alt}"
        status[rid] = ("resolved", ";".join(notes), vid)
        sites.append((contig_rank[chrom], pos, chrom, f"{ref},{alt}"))
        scores.append({"ID": vid, "EFFECT_ALLELE": ea, "OTHER_ALLELE": oa,
                       "WEIGHT": info.loc[rid, "Log_Odds_Ratio"],
                       "EAF": info.loc[rid, "Effect_Allele_Frequency"],
                       "REF": ref, "ALT": alt, "row_id": rid})

    scores_df = pd.DataFrame(scores)
    dup = scores_df["ID"].duplicated(keep=False) if not scores_df.empty else pd.Series(dtype=bool)
    if dup.any():
        sys.exit(f"Duplicate variant IDs after normalisation:\n{scores_df[dup]}")

    sites.sort()
    with open(os.path.join(args.outdir, "sites.tsv"), "w") as fh:
        for _, pos, chrom, alleles in sites:
            fh.write(f"{chrom}\t{pos}\t{alleles}\n")
    scores_df.to_csv(os.path.join(args.outdir, "score_plink2.tsv"), sep="\t", index=False)

    st = df[["row_id"] + CANRISK_COLS].copy()
    st["status"] = st["row_id"].map(lambda r: status[r][0])
    st["detail"] = st["row_id"].map(lambda r: status[r][1])
    st["final_id"] = st["row_id"].map(lambda r: status[r][2] if len(status[r]) > 2 else "")
    st["pgs_check"] = st["row_id"].map(lambda r: pgs_note.get(r, ""))
    st.to_csv(os.path.join(args.outdir, "site_status.tsv"), sep="\t", index=False)
    with open(os.path.join(args.outdir, "prs_info.tsv"), "w") as fh:
        fh.write(f"source\t{os.path.basename(args.canrisk)}\nalpha\t{alpha}\nn_variants\t{len(df)}\n")

    print(f"Source                      : {args.canrisk} (alpha = {alpha})")
    print(f"Total variants              : {len(df)}")
    print(f"Resolved                    : {(st['status'] == 'resolved').sum()}")
    print(f"Unresolved                  : {(st['status'] == 'unresolved').sum()}")
    print("Detail counts:")
    print(st["detail"].str.replace(r":[^;]*", "", regex=True).value_counts().to_string())
    if pgs_pos is not None:
        print("PGS Catalog position cross-check:")
        print(st["pgs_check"].str.split(":").str[0].value_counts().to_string())


if __name__ == "__main__":
    main()
