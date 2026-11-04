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
The only required inputs are the alignments (indexed BAM or CRAM) and the reference FASTA they
were aligned to (with `.fai`). The pipeline first checks that every panel REF allele matches this
FASTA and stops if any does not.

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
