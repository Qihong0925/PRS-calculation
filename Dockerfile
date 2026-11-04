# Container for the PRS-from-BAM pipeline. Tools and versions come from env/nextflow_env.yml.
FROM mambaorg/micromamba:1.5.10-bookworm-slim

LABEL org.opencontainers.image.title="prs-from-bam" \
      org.opencontainers.image.description="bcftools/samtools/Python environment for the PRS-from-BAM Nextflow pipeline"

COPY --chown=$MAMBA_USER:$MAMBA_USER env/nextflow_env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean -a -y && \
    rm /tmp/env.yml

ENV PATH="/opt/conda/bin:$PATH" \
    PYTHONNOUSERSITE=1 \
    MPLBACKEND=Agg
