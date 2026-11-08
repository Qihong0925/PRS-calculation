# PRS from BAM (GRCh38)

Nextflow pipeline that computes breast cancer polygenic risk scores (CanRisk BCAC 313 and
BRIDGES 306) from germline WGS BAM/CRAM files. Panel sites are force-genotyped with bcftools,
missing sites are imputed at 2 x EAF, and scores are standardised to the CanRisk population scale.

## Quick start
```bash
git clone https://github.com/Qihong0925/PRS-calculation.git
cd PRS-calculation

# BAM/CRAM glob (sample name = file name without .bam/.cram, or use --strip_suffix)
nextflow run main.nf -profile slurm,singularity \
    --input '/path/to/bams/*.bam' \
    --fasta /path/to/GRCh38.fa \
    --batch MY_BATCH --outdir results

# or a samplesheet (.csv/.tsv) with columns: sample, bam[, index]
nextflow run main.nf -profile slurm,singularity --input samples.tsv --fasta /path/to/GRCh38.fa
```
The only required inputs are **germline** alignments (indexed BAM or CRAM) and the **GRCh38**
reference FASTA they were aligned to (chr-prefixed contigs, with `.fai`; ALT contigs are fine).
The pipeline stops with a plain message if the FASTA is GRCh37/hg19 or lacks the `chr` prefix, and
checks that every panel REF allele matches the FASTA. Run `nextflow run main.nf --help` for all options.

Test on a tiny synthetic data set (about 1 min):
```bash
nextflow run main.nf -profile test,singularity   # or test,docker / test,conda
```

## Options
| Parameter | Default | Description |
|---|---|---|
| `--input` | required | BAM/CRAM glob, or samplesheet with columns `sample`, `bam` or `cram` (optional `index`) |
| `--fasta` | required | Alignment reference FASTA (with `.fai`) |
| `--outdir` | `results` | Output directory |
| `--batch` | `cohort` | Label used in output file names |
| `--panels` | `BCAC_313,BRIDGES_306` | Panels under `--panel_dir`; first = primary, second = secondary |
| `--strip_prefix` / `--strip_suffix` | none | Trim file names to sample names (glob input) |
| `--min_mapq` / `--min_baseq` | 20 / 20 | Read and base quality filters |
| `--force_expected` | none | TSV (column `site` CHROM:POS or `ID`) of sites not trusted in this batch; scored at 2 x EAF |
| `--tumour_sheet` | none | Tumour mode instead of `--input` (see below) |
| `--queue` | cluster default | SLURM partition for all jobs (with `-profile slurm`) |
| `--genotype_max_forks` | 0 (no limit) | Max concurrent genotyping jobs, e.g. 1 for files on tiered/archived storage |
| `--genotype_retries` / `--retry_wait_sec` | 2 / 0 | Retries per sample and wait between them |

Profiles: `slurm` (SLURM executor), `singularity` / `docker` (pipeline container, recommended),
`conda` (tools from `env/nextflow_env.yml`), `env_prs` (use `env/env_prs.sh` + `config/local.sh`:
environment modules and an existing conda env), `test`. Site-specific executor settings
(queue, account) can go in a `*.local.config` passed with `-c`.

## Software environment
All tools are defined in `env/nextflow_env.yml` (bcftools 1.19, htslib 1.19.1, samtools 1.19.2,
Python 3.11 with numpy 1.26.4, pandas 2.2.3, scipy 1.17.1). The same file is used for `-profile conda`
and to build the container (`Dockerfile`). GitHub Actions (`.github/workflows/container.yml`)
builds the image on every push, publishes it to `ghcr.io/qihong0925/prs-from-bam`, and runs the test
profile inside it.
```bash
nextflow run main.nf -profile slurm,singularity --container ghcr.io/qihong0925/prs-from-bam:latest \
    --input '/path/to/bams/*.bam' --fasta /path/to/GRCh38.fa
```

## Tumour mode: tumour vs germline PRS
Compares each tumour with its matched normal to test whether tumours preferentially retain
risk or protective alleles at PRS sites (loss of heterozygosity / allelic imbalance).
```bash
nextflow run main.nf -profile slurm,singularity --tumour_sheet tumours.tsv \
    --fasta /path/to/GRCh38.fa --batch MY_BATCH --outdir results_tumour
```
`tumours.tsv` (tab- or comma-separated) has the columns `tumour`, `normal`, `tumour_path`,
`normal_path` and `purity` (fraction, e.g. 0.45). A normal can be listed with several tumours.
Relative paths are resolved against the sheet's folder. Tumours with purity below `--min_purity`
(0.15) are listed in `qc/tumour_excluded.tsv` and not analysed. Copy number is not used.

Method, at each germline-heterozygous PRS site with tumour depth >= `--tumour_min_dp` (20), using
the tumour effect-allele fraction v and purity p:

| Tumour effect-allele VAF | Call | Tumour-cell effect dosage |
|---|---|---|
| v >= 0.70 | effect allele retained (any purity) | 2 |
| 0.55 <= v < 0.70 | retained if v >= 0.5 + 0.8 x p/2, otherwise no change | 2 or 1 |
| 0.45 < v < 0.55 | no change (any purity) | 1 |
| 0.30 < v <= 0.45 | lost if v <= 0.5 - 0.8 x p/2, otherwise no change | 0 or 1 |
| v <= 0.30 | effect allele lost (any purity) | 0 |

p/2 is the largest VAF shift complete LOH can produce at purity p; 0.8 (`--loh_purity_fraction`)
requires 80% of it. Homozygous, missing, forced and low-depth sites keep the germline dosage.
The tumour-cell PRS uses these dosages; the germline PRS is the normal's score. The risk allele is
the effect allele when its weight is positive, otherwise the other allele.

Outputs in `tumour/`:
- `<batch>.tumour_prs.tsv`: per tumour and panel, germline / tumour-cell / directly called scores
  (raw, z, percentile), their differences, and counts of risk- vs protective-allele retention;
- `<batch>.tumour_sites.tsv`: per tumour, panel and site, the genotypes, depths, VAFs and calls;
- `<batch>.tumour_summary.tsv`: cohort tests per panel (paired Wilcoxon and sign test on
  tumour-cell minus germline z; binomial test of risk vs protective retention; Wilcoxon test of the
  tumour-minus-normal risk-allele VAF shift), for all tumours and for one tumour per patient
  (highest purity);
- `<batch>.tumour_vs_germline.png`.

Test: `nextflow run main.nf -profile test_tumour,singularity` (simulated tumours with known allele
loss; `tests/data/expected_tumour_events.tsv` holds the truth).

## Interpreting the scores
- **Use germline samples only.** The pipeline cannot tell tumour from normal: a tumour BAM produces
  a score without any warning, but somatic copy-number changes and LOH distort the genotypes.
- **Primary score**: `<panel>_z` / `<panel>_percentile` in `final_scores.tsv` (BCAC 313 is the
  primary panel). Missing sites are filled with the expected dosage 2 x EAF; `n_missing` counts them.
- **Standardisation**: z and percentiles use the CanRisk population mean and SD
  (`resources/canrisk_standardisation.tsv`), which were derived in women of **European ancestry**.
  For other ancestries the raw score ranks samples, but the percentile is not calibrated:
  check ancestry before reporting percentiles.
- **Sensitivity columns**: `_bam` scores use the calls exactly as genotyped; `_flagged_as_expected`
  sets validation-flagged sites (`panels/*/validation/site_validation.flagged.tsv`) to 2 x EAF;
  the secondary panel (BRIDGES 306) gives an independent estimate. Large differences for a sample
  point to genotyping problems at a few sites (see `genotypes/<sample>.gt.tsv`).
- A PRS is a research measure here; it is not a clinical risk estimate on its own.

## Outputs
| Path | Contents |
|---|---|
| `final/<batch>.final_scores.tsv` | Per-sample raw score, z and percentile for each panel (corrected and direct-from-BAM) |
| `final/<batch>.sensitivity.tsv` | Primary vs flagged-as-expected vs secondary panel, and corrected vs BAM |
| `final/<batch>.cohort_summary.tsv`, `percentile_bands.tsv`, `*.png/pdf` | Cohort distribution vs population |
| `scores/<panel>.scores.tsv`, `.site_summary.tsv`, `.missing_sites.tsv` | Per-panel scores and per-site call rates / EAF |
| `scores/cohort.prs_sites.vcf.gz` | Cohort VCF of panel-constrained calls |
| `genotypes/<sample>.*` | Constrained and unconstrained VCFs, per-site GT/DP/AD/GQ, no-coverage sites, log |
| `qc/` | Reference check and the sample manifest used |

## Pipeline steps
`00` reference check -> `03a` union of panel sites -> `03` per-sample genotyping (in parallel) ->
`04` cohort merge and scoring -> `05` sensitivity analysis -> `05b` cohort summary.
Each step wraps the scripts in `bin/` and `py/`. `bin/03_genotype_array.sh` and
`bin/03_genotype_sequential.sh` are the earlier stand-alone SLURM runners and are not used by the pipeline.

## Panels (steps 01-02, run once)
The frozen panels in `panels/` were built with `build_panel.sh`. To rebuild:
```bash
conda env create -n PRS -f env/prs_env.lock.yml
cp config/local.example.sh config/local.sh   # set CanRisk, chain, FASTA, gnomAD paths
PRS_NAME=BCAC_313 bash build_panel.sh
```
Built: 2026-10-27 by Qihong.

## Sources
- Weights: CanRisk knowledge base, "What variants are used in the PRS?"
  (BCAC_313_PRS.prs, BRIDGES_306_PRS.prs; GRCh37; downloaded 2026-10-27)
- Liftover: UCSC hg19ToHg38.over.chain.gz (pyliftover)
- Reference: GRCh38.p13 full assembly including ALT contigs (used for REF-allele checks)
- Position cross-check: PGS Catalog PGS000004 hmPOS_GRCh38 (HmPOS_date 2022-07-29)
- Validation: gnomAD v4.1 genomes (AF_nfe)

## Method (steps 01-02)
CanRisk REF/effect alleles are taken as truth. GRCh37 -> GRCh38 by chain liftover;
SNVs on '-' strand reverse-complemented; '-' strand indels require a documented manual
correction (resources/manual_corrections.tsv). REF checked against GRCh38 FASTA; sites where the
GRCh38 reference carries the effect allele are swapped with effect allele preserved.
All variants left-aligned with bcftools norm. No variant is dropped; site_status.tsv lists all.

## BCAC_313 summary
313/313 resolved; 1 manual correction (rs143384623, 1q21 inversion); 2 normalised indels;
1 GRCh38 REF = effect allele (chr9:133271182). gnomAD: 310 exact, 3 absent;
15 flagged (validation/site_validation.flagged.tsv): multi-allelic STRs, known NGS-problem
loci (several removed/replaced in BRIDGES 306), and poorly callable regions. All coordinates of
flagged SNVs confirmed via rsID GRCh38 mapping (Ensembl). Flagged sites are kept in the primary
score and handled by sensitivity analysis.
