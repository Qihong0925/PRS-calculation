#!/usr/bin/env nextflow
/*
 * PRS-from-BAM: force-genotype frozen PRS panel sites in germline BAM/CRAM files and score them.
 *
 *   00 CHECK_REFERENCE  panel REF alleles must match the alignment reference FASTA
 *   03a UNION_SITES     union of panel sites (each sample is genotyped once for all panels)
 *   03 GENOTYPE         per-sample constrained calls at panel sites (bcftools mpileup/call)
 *   04 MERGE_SCORE      cohort VCF + per-panel scores (observed / imputed / flagged / BAM)
 *   05 SENSITIVITY      primary vs secondary panel and corrected vs BAM scores; final table
 *   05b SUMMARY         cohort summary tables and figures
 *
 * Panels are built once with build_panel.sh (steps 01-02) and shipped in panels/.
 */

def sampleName(f) {
    def name = f.name
    if (params.strip_prefix && name.startsWith(params.strip_prefix)) {
        name = name.substring(params.strip_prefix.length())
    }
    if (params.strip_suffix && name.endsWith(params.strip_suffix)) {
        name = name.substring(0, name.length() - params.strip_suffix.length())
    } else {
        name = name.replaceAll(/\.(bam|cram)$/, '')
    }
    return name
}

def findIndex(f) {
    def base = f.toString().replaceAll(/\.(bam|cram)$/, '')
    def candidates = ["${f}.bai", "${f}.crai", "${base}.bai", "${base}.crai", "${f}.csi"]
    def found = candidates.collect { c -> file(c) }.find { idx -> idx.exists() }
    if (!found) {
        error "No index (.bai/.crai/.csi) found for ${f}"
    }
    return found
}

def helpMessage() {
    return """
    PRS-from-BAM: breast cancer polygenic risk scores (CanRisk BCAC 313 / BRIDGES 306) from germline WGS

    Usage:
      nextflow run Qihong0925/PRS-calculation -profile <slurm,>singularity \\
          --input '/path/to/bams/*.bam' --fasta /path/to/GRCh38.fa --batch MY_BATCH --outdir results

    Required:
      --input            BAM/CRAM glob (quote it), or samplesheet .csv/.tsv with columns
                         'sample' and 'bam' (or 'cram'), optional 'index'. GERMLINE samples only.
      --fasta            GRCh38 FASTA the alignments were made against (chr-prefixed contigs, with .fai)

    Optional:
      --outdir           Output directory [${params.outdir}]
      --batch            Label used in output file names [${params.batch}]
      --panels           Comma-separated panels in --panel_dir; first = primary [${params.panels}]
      --strip_prefix     Remove this prefix from file names to get sample names (glob input)
      --strip_suffix     Remove this suffix instead of .bam/.cram (e.g. .dedup.bqsr.bam)
      --min_mapq         Minimum read mapping quality [${params.min_mapq}]
      --min_baseq        Minimum base quality [${params.min_baseq}]
      --force_expected   TSV (column 'site' CHROM:POS or 'ID') of sites scored at 2 x EAF
      --queue            SLURM partition for all jobs (with -profile slurm)
      --genotype_max_forks  Max concurrent genotyping jobs; e.g. 3 for archived storage (0 = no limit)
      --genotype_retries    Retries per sample [${params.genotype_retries}]
      --retry_wait_sec      Seconds to wait before a retry [${params.retry_wait_sec}]

    Profiles (-profile, comma-separated):
      singularity | docker | conda   software (container recommended)
      slurm                          submit jobs to SLURM
      test                           tiny synthetic data set (about 1 min)
    """.stripIndent()
}

def checkReference(fai) {
    def len = [:]
    fai.readLines().each { line ->
        def f = line.tokenize('\t')
        len[f[0]] = f[1] as long
    }
    if (!len.containsKey('chr1')) {
        if (len.containsKey('1')) {
            error "Reference ${params.fasta} has contigs named '1', '2', ... The panels use GRCh38 with 'chr' prefixes (chr1, chr2, ...). Use the chr-prefixed GRCh38 FASTA your BAMs were aligned to."
        }
        error "Reference ${params.fasta} has no contig 'chr1'. A GRCh38 reference with chr-prefixed contigs is required."
    }
    if (len['chr1'] == 249250621L) {
        error "Reference ${params.fasta} looks like GRCh37/hg19 (chr1 length 249,250,621). The panels are GRCh38: re-align to GRCh38 or use the GRCh38 FASTA."
    }
    if (len['chr1'] != 248956422L) {
        error "Reference ${params.fasta}: chr1 length ${len['chr1']} is not GRCh38 (248,956,422)."
    }
    def n_alt = len.keySet().count { k -> k.endsWith('_alt') }
    log.info "Reference: GRCh38, ${len.size()} contigs (${n_alt} ALT)"
}

process CHECK_REFERENCE {
    tag "${fasta.name}"
    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    tuple path(fasta), path(fai)
    path panels

    output:
    path 'reference_check.txt'

    script:
    """
    bash ${projectDir}/bin/00_check_reference.sh ${fasta} ${panels} 2>&1 | tee reference_check.txt
    """

    stub:
    """
    echo "stub" > reference_check.txt
    """
}

process UNION_SITES {
    publishDir "${params.outdir}/sites", mode: params.publish_mode

    input:
    path panels
    path reference_check

    output:
    tuple path('union_sites.tsv.gz'), path('union_sites.tsv.gz.tbi'), path('union_sites.regions.tsv.gz'), path('union_sites.regions.tsv.gz.tbi')

    script:
    """
    bash ${projectDir}/bin/03a_union_sites.sh union_sites.tsv.gz ${panels}
    """

    stub:
    """
    touch union_sites.tsv.gz union_sites.tsv.gz.tbi union_sites.regions.tsv.gz union_sites.regions.tsv.gz.tbi
    """
}

process GENOTYPE {
    tag "${sample}"
    publishDir "${params.outdir}/genotypes", mode: params.publish_mode

    input:
    tuple val(sample), path(bam), path(bam_index)
    tuple path(fasta), path(fai)
    tuple path(sites), path(sites_tbi), path(regions), path(regions_tbi)

    output:
    tuple val(sample), path("${sample}.vcf.gz"), path("${sample}.vcf.gz.tbi"), emit: vcf
    path "${sample}.free.vcf.gz*", emit: free
    path "${sample}.{gt,nocov}.tsv", emit: tables
    path "${sample}.genotype.log", emit: log

    script:
    def ext = bam_index.name.endsWith('.crai') ? 'crai' : (bam_index.name.endsWith('.csi') ? 'csi' : 'bai')
    """
    # htslib looks for <alignment>.<ext>; link the staged index under that name if needed
    [[ -e "${bam}.${ext}" ]] || ln -s "${bam_index}" "${bam}.${ext}"
    # Read the index first: forces recall of files on tiered/archived storage
    cat "${bam_index}" > /dev/null
    bash ${projectDir}/bin/03_genotype_sample.sh "${sample}" "${bam}" "${fasta}" "${sites}" . \\
        ${params.min_mapq} ${params.min_baseq} 2>&1 | tee "${sample}.genotype.log"
    """

    stub:
    """
    touch ${sample}.vcf.gz ${sample}.vcf.gz.tbi ${sample}.free.vcf.gz ${sample}.free.vcf.gz.tbi
    touch ${sample}.gt.tsv ${sample}.nocov.tsv ${sample}.genotype.log
    """
}

process MERGE_SCORE {
    publishDir "${params.outdir}", mode: params.publish_mode

    input:
    val samples
    path vcfs, stageAs: 'gt/*'
    path panels
    path force_expected

    output:
    path 'scores', emit: scores

    script:
    def fe = force_expected.name != 'NO_FILE' ? "FORCE_EXPECTED=${force_expected}" : ''
    """
    printf "%s\\n" sample ${samples.join(' ')} > samples.tsv
    ${fe} bash ${projectDir}/bin/04_merge_score.sh ${projectDir} samples.tsv gt scores ${panels}
    """

    stub:
    """
    mkdir -p scores
    for p in ${panels}; do touch scores/\$(basename \$p).scores.tsv; done
    touch scores/cohort.prs_sites.vcf.gz
    """
}

process SENSITIVITY {
    publishDir "${params.outdir}/final", mode: params.publish_mode

    input:
    path scores
    tuple val(primary), val(secondary)

    output:
    path "${params.batch}.final_scores.tsv", emit: final_scores
    path "${params.batch}.sensitivity.tsv", emit: sensitivity

    script:
    """
    python ${projectDir}/py/sensitivity.py \\
        --score-dir ${scores} \\
        --primary ${primary} \\
        --secondary ${secondary} \\
        --standardisation ${projectDir}/resources/canrisk_standardisation.tsv \\
        --out-prefix ${params.batch}
    """

    stub:
    """
    touch ${params.batch}.final_scores.tsv ${params.batch}.sensitivity.tsv
    """
}

process SUMMARY {
    publishDir "${params.outdir}/final", mode: params.publish_mode

    input:
    path final_scores
    tuple val(primary), val(secondary)

    output:
    path "${params.batch}.*"

    script:
    """
    python ${projectDir}/py/cohort_summary.py \\
        --final-scores ${final_scores} \\
        --primary ${primary} \\
        --secondary ${secondary} \\
        --batch ${params.batch} \\
        --out-prefix ${params.batch}
    python ${projectDir}/py/plot_prs_distribution.py ${final_scores} ${params.batch} ${params.batch}.kde ${primary} ${secondary}
    """

    stub:
    """
    touch ${params.batch}.cohort_summary.tsv ${params.batch}.kde.prs_distribution.png
    """
}

workflow {
    main:
    if (params.help) {
        log.info helpMessage()
        return
    }
    if (!params.input) {
        error "Missing --input (BAM/CRAM glob, or samplesheet .csv/.tsv with columns 'sample' and 'bam'). Run with --help for usage."
    }
    if (!params.fasta) {
        error "Missing --fasta (GRCh38 FASTA the BAMs were aligned to; needs a .fai index). Run with --help for usage."
    }

    // Reference
    def fasta = file(params.fasta, checkIfExists: true)
    def fai = file("${params.fasta}.fai")
    if (!fai.exists()) {
        error "Missing FASTA index ${fai}; run: samtools faidx ${params.fasta}"
    }
    if (params.check_genome) {
        checkReference(fai)
    }
    ch_ref = channel.value([fasta, fai])

    // Panels (first = primary, second = secondary)
    def panel_names = params.panels.toString().tokenize(',')*.trim()
    def panel_dirs = panel_names.collect { p -> file("${params.panel_dir}/${p}", checkIfExists: true) }
    ch_panels = channel.value(panel_dirs)

    // Alignments: samplesheet (columns sample, bam[, index]) or a glob of BAM/CRAM files
    def input = params.input.toString()
    if (input.endsWith('.csv') || input.endsWith('.tsv')) {
        ch_bams = channel.fromPath(input, checkIfExists: true)
            .splitCsv(header: true, sep: input.endsWith('.tsv') ? '\t' : ',')
            .map { row ->
                def aln = row.bam ?: row.cram
                if (!row.sample || !aln) {
                    error "Samplesheet ${input} needs columns 'sample' and 'bam' (or 'cram')"
                }
                def bam = file(aln, checkIfExists: true)
                def idx = row.index ? file(row.index, checkIfExists: true) : findIndex(bam)
                [row.sample, bam, idx]
            }
    } else {
        ch_bams = channel.fromPath(input, checkIfExists: true)
            .map { bam -> [sampleName(bam), bam, findIndex(bam)] }
    }

    // Fail on duplicate sample names; record the sample manifest
    ch_bams = ch_bams
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { rows ->
            def names = rows.collect { r -> r[0] }
            def dups = names.findAll { n -> names.count(n) > 1 }.unique()
            if (dups) {
                error "Duplicate sample names: ${dups.join(', ')}"
            }
            log.info "PRS: ${rows.size()} samples; panels ${panel_names.join(', ')}"
            rows
        }
        .flatMap { rows -> rows }

    ch_bams
        .map { s, b, i -> "${s}\t${b}\t${i}\n" }
        .collectFile(name: 'samples.tsv', seed: "sample\tbam\tindex\n", sort: true, storeDir: "${params.outdir}/qc")

    // 00 / 03a / 03
    CHECK_REFERENCE(ch_ref, ch_panels)
    UNION_SITES(ch_panels, CHECK_REFERENCE.out)
    GENOTYPE(ch_bams, ch_ref, UNION_SITES.out)

    // 04
    ch_gt = GENOTYPE.out.vcf
        .toSortedList { a, b -> a[0] <=> b[0] }
        .multiMap { rows ->
            names: rows.collect { r -> r[0] }
            files: rows.collect { r -> [r[1], r[2]] }.flatten()
        }
    def fe = params.force_expected ? file(params.force_expected, checkIfExists: true) : file("${projectDir}/assets/NO_FILE")
    MERGE_SCORE(ch_gt.names, ch_gt.files, ch_panels, fe)

    // 05 / 05b (need a primary and a secondary panel)
    if (panel_names.size() >= 2) {
        ch_pair = channel.value([panel_names[0], panel_names[1]])
        SENSITIVITY(MERGE_SCORE.out.scores, ch_pair)
        SUMMARY(SENSITIVITY.out.final_scores, ch_pair)
    } else {
        log.warn "Only one panel given: skipping sensitivity analysis and cohort summary"
    }
}
