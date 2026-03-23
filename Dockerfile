# Base image — mambaforge provides conda + mamba pre-installed on a
# minimal Ubuntu base. Pinned to a specific version for reproducibility.
FROM condaforge/mambaforge:23.3.1-1

# Metadata
LABEL maintainer="Scott Lewis"
LABEL description="exRNA-NF: sRNA diversity pipeline"
LABEL version="1.0.0"

# Set working directory
WORKDIR /pipeline

# Copy the conda environment dependency specification
COPY environment.yml .

# Install all conda/mamba dependencies from environment.yml in one layer.
# --no-capture-output streams install progress to the build log.
# clean removes cached packages after install to reduce image size.
RUN mamba env create -f environment.yml && \
    mamba clean --all --yes

# Make conda environment activation automatic for all subsequent RUN,
# CMD, and ENTRYPOINT commands by prepending the env bin to PATH.
# Replace 'ex-srna-nf' with whatever name is set in environment.yml.
ENV PATH="/opt/conda/envs/exRNA/bin:$PATH"

# Build pipeline source files into the image
COPY main_sRNA.nf .
COPY nextflow.config .
COPY modules/ modules/
COPY bin/ bin/

# Ensure bin/ scripts are executable
RUN chmod +x bin/*.py

# Verify tool availability at build time
RUN nextflow -version && \
    bowtie2 --version && \
    samtools --version && \
    python3 -c "import pysam; print('pysam ok')"

# Default command — print usage when container is run with no arguments
CMD ["nextflow", "run", "main_sRNA.nf", "--help"]
