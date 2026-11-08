#!/usr/bin/env python3
"""Generate a tiny synthetic test set for the Nextflow pipeline (-profile test).

Creates under tests/data/:
  ref.fa(.fai)                 one 20 kb random contig 'chrT'
  panels/BCAC_313, BRIDGES_306 mini panels (same file layout as real panels; names match
                               resources/canrisk_standardisation.tsv so steps 05/05b run)
  bams/                        S1 (.bam + .bam.bai), S2 (.bam + .bai), S3 (.cram + .crai)
  expected_genotypes.tsv       simulated truth (alt-allele count per sample and site)
  bams/S1T, S2T, S3T           tumours of S1-S3 (purity 0.6, 0.3, 0.1) with known allele loss
  tumour_sheet.tsv             tumour sheet for -profile test_tumour (paths relative to it)
  expected_tumour_events.tsv   simulated tumour-cell ALT fraction at each germline-heterozygous site
One panel site has no reads in any sample (exercises the no-coverage path). Tumour data are generated
after the germline data, so the germline files do not depend on them.
Requires samtools, bgzip and tabix on PATH.
"""
import os
import random
import subprocess

random.seed(1)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")
L = 20000
READ_LEN, DEPTH = 100, 30
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def run(cmd):
    subprocess.run(cmd, shell=True, check=True)


os.makedirs(os.path.join(OUT, "bams"), exist_ok=True)
seq = "".join(random.choice("ACGT") for _ in range(L))
with open(os.path.join(OUT, "ref.fa"), "w") as fh:
    fh.write(">chrT\n" + "\n".join(seq[i:i + 60] for i in range(0, L, 60)) + "\n")
run(f"samtools faidx {OUT}/ref.fa")

positions = [1000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 18000]   # 1-based
nocov = 18000
sites = []
for i, pos in enumerate(positions):
    ref = seq[pos - 1]
    alt = random.choice([b for b in "ACGT" if b != ref])
    sites.append(dict(pos=pos, ref=ref, alt=alt, eaf=round(random.uniform(0.1, 0.6), 3),
                      w=round(random.uniform(-0.1, 0.1), 4), ea_is_alt=(i % 3 != 0)))

panels = {"BCAC_313": sites, "BRIDGES_306": sites[1:]}
for name, ps in panels.items():
    d = os.path.join(OUT, "panels", name)
    os.makedirs(os.path.join(d, "validation"), exist_ok=True)
    with open(os.path.join(d, "sites.tsv"), "w") as fh:
        for s in ps:
            fh.write(f"chrT\t{s['pos']}\t{s['ref']},{s['alt']}\n")
    run(f"bgzip -f {d}/sites.tsv && tabix -f -s1 -b2 -e2 {d}/sites.tsv.gz")
    with open(os.path.join(d, "score_plink2.tsv"), "w") as fh:
        fh.write("ID\tEFFECT_ALLELE\tOTHER_ALLELE\tWEIGHT\tEAF\tREF\tALT\trow_id\n")
        for i, s in enumerate(ps):
            ea, oa = (s["alt"], s["ref"]) if s["ea_is_alt"] else (s["ref"], s["alt"])
            fh.write(f"chrT:{s['pos']}:{s['ref']}:{s['alt']}\t{ea}\t{oa}\t{s['w']}\t{s['eaf']}\t{s['ref']}\t{s['alt']}\tv{i:04d}\n")
    with open(os.path.join(d, "validation", "site_validation.flagged.tsv"), "w") as fh:
        s = ps[1]
        fh.write(f"row_id\tID\nv0001\tchrT:{s['pos']}:{s['ref']}:{s['alt']}\n")

samples = {"S1": ("bam", ".bam.bai"), "S2": ("bam", ".bai"), "S3": ("cram", ".crai")}
truth = []
for sname, (fmt, idx_ext) in samples.items():
    sam = os.path.join(OUT, f"{sname}.sam")
    with open(sam, "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:{L}\n@RG\tID:{sname}\tSM:{sname}\n")
        n = 0
        for s in sites:
            g = random.choice([0, 1, 2])
            truth.append((sname, s["pos"], g))
            if s["pos"] == nocov:
                continue
            for k in range(DEPTH):
                start = s["pos"] - random.randint(10, READ_LEN - 10)     # 1-based leftmost
                r = list(seq[start - 1:start - 1 + READ_LEN])
                carries_alt = (g == 2) or (g == 1 and k % 2 == 0)
                if carries_alt:
                    r[s["pos"] - start] = s["alt"]
                n += 1
                fh.write(f"r{n}\t0\tchrT\t{start}\t60\t{READ_LEN}M\t*\t0\t0\t{''.join(r)}\t{'I' * READ_LEN}\tRG:Z:{sname}\n")
    base = os.path.join(OUT, "bams", sname)
    if fmt == "bam":
        run(f"samtools sort -o {base}.bam {sam}")
        run(f"samtools index -o {base}{idx_ext} {base}.bam")
    else:
        run(f"samtools sort -O cram --reference {OUT}/ref.fa -o {base}.cram {sam}")
        run(f"samtools index {base}.cram")
    os.remove(sam)

with open(os.path.join(OUT, "expected_genotypes.tsv"), "w") as fh:
    fh.write("sample\tpos\talt_count\n")
    for t in truth:
        fh.write("\t".join(map(str, t)) + "\n")
# ---- Tumours: tumour-cell ALT fraction at germline-heterozygous sites is set per site ----
# Reads are drawn deterministically: ALT reads = round(depth x VAF),
# VAF = purity x f_tumour + (1 - purity) x 0.5
T_DEPTH, LOW_DEPTH = 80, 12
tumours = {"S1T": ("S1", 0.6), "S2T": ("S2", 0.3), "S3T": ("S3", 0.1)}
events = ["keep_alt", "keep_ref", "partial", "none", "low_depth"]
f_tumour = {"keep_alt": 1.0, "keep_ref": 0.0, "partial": 0.75, "none": 0.5, "low_depth": 1.0}
geno = {(sm, pos): g for sm, pos, g in truth}
tevents, k_event = [], 0
for tname, (normal, purity) in tumours.items():
    sam = os.path.join(OUT, f"{tname}.sam")
    with open(sam, "w") as fh:
        fh.write(f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:{L}\n@RG\tID:{tname}\tSM:{tname}\n")
        n = 0
        for s in sites:
            g = geno[(normal, s["pos"])]
            if s["pos"] == nocov:
                continue
            depth = T_DEPTH
            if g == 1:
                ev = events[k_event % len(events)]
                k_event += 1
                if ev == "low_depth":
                    depth = LOW_DEPTH
                vaf = purity * f_tumour[ev] + (1 - purity) * 0.5
                tevents.append((tname, normal, purity, s["pos"], ev, f_tumour[ev], depth, round(depth * vaf)))
                n_alt = round(depth * vaf)
            else:
                n_alt = depth if g == 2 else 0
            for k in range(depth):
                start = s["pos"] - random.randint(10, READ_LEN - 10)
                r = list(seq[start - 1:start - 1 + READ_LEN])
                if k < n_alt:
                    r[s["pos"] - start] = s["alt"]
                n += 1
                fh.write(f"t{n}\t0\tchrT\t{start}\t60\t{READ_LEN}M\t*\t0\t0\t{''.join(r)}\t{'I' * READ_LEN}\tRG:Z:{tname}\n")
    base = os.path.join(OUT, "bams", tname)
    run(f"samtools sort -o {base}.bam {sam}")
    run(f"samtools index {base}.bam")
    os.remove(sam)

with open(os.path.join(OUT, "expected_tumour_events.tsv"), "w") as fh:
    fh.write("tumour\tnormal\tpurity\tpos\tevent\ttumour_cell_alt_fraction\tdepth\talt_reads\n")
    for t in tevents:
        fh.write("\t".join(map(str, t)) + "\n")
normal_path = {"S1": "bams/S1.bam", "S2": "bams/S2.bam", "S3": "bams/S3.cram"}
with open(os.path.join(OUT, "tumour_sheet.tsv"), "w") as fh:
    fh.write("tumour\tnormal\ttumour_path\tnormal_path\tpurity\n")
    for tname, (normal, purity) in tumours.items():
        fh.write(f"{tname}\t{normal}\tbams/{tname}.bam\t{normal_path[normal]}\t{purity}\n")
print("Test data written to", OUT)
