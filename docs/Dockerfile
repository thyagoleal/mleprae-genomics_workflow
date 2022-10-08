
=======================================================================
[95m>>> Processing 3 samples in SE format:
[0m
[92mMM20743-100k-N-reads
MM114-100k
MM23897-100k-N-reads
[0m
[94m>>> Processing 3 samples in PE format:
[0m
[92mSRR6241727-100k
SRR6241756-100k
MM27793-100k[0m
=======================================================================
FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="043181009c60f4699ac77875ea60692fd2a77d0920bb75faea0a425b18bccac1"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/basic.yml
#   prefix: /conda-envs/3bf9c50269bab91f237f6667f83ab34a
#   name: python basics
#   channels:
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - python=3.10.5
#     - pandas=1.4.3
#     - snakemake=7.5.0
#     - snakemake-minimal=7.5.0
#     - snakemake-wrapper-utils=0.5.0
RUN mkdir -p /conda-envs/3bf9c50269bab91f237f6667f83ab34a
COPY envs/basic.yml /conda-envs/3bf9c50269bab91f237f6667f83ab34a/environment.yaml

# Conda environment:
#   source: envs/mapping_tools.yml
#   prefix: /conda-envs/e5d860111254c639477a4e8398f307bf
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bowtie2=2.4.5
#     - samtools=1.15
#     - bcftools=1.15
#     - tabixpp=1.1.0
#     - freebayes=1.3.6
RUN mkdir -p /conda-envs/e5d860111254c639477a4e8398f307bf
COPY envs/mapping_tools.yml /conda-envs/e5d860111254c639477a4e8398f307bf/environment.yaml

# Conda environment:
#   source: envs/picard.yml
#   prefix: /conda-envs/f748ab4e3674d4b980a19edf89a39484
#   channels:
#    - bioconda
#   dependencies:
#    - picard=2.27.4
RUN mkdir -p /conda-envs/f748ab4e3674d4b980a19edf89a39484
COPY envs/picard.yml /conda-envs/f748ab4e3674d4b980a19edf89a39484/environment.yaml

# Conda environment:
#   source: envs/qc_trimming.yml
#   prefix: /conda-envs/abfc4906b35726017b76a70cae2c256f
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - fastqc=0.11.9
#     - trimmomatic=0.39
RUN mkdir -p /conda-envs/abfc4906b35726017b76a70cae2c256f
COPY envs/qc_trimming.yml /conda-envs/abfc4906b35726017b76a70cae2c256f/environment.yaml

# Conda environment:
#   source: envs/qualimap.yml
#   prefix: /conda-envs/f03630632771bbe472a0c60a72c3b52b
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - qualimap=2.2.2d
RUN mkdir -p /conda-envs/f03630632771bbe472a0c60a72c3b52b
COPY envs/qualimap.yml /conda-envs/f03630632771bbe472a0c60a72c3b52b/environment.yaml

# Conda environment:
#   source: envs/seqprep.yml
#   prefix: /conda-envs/1f131affa68a486db397c65b7a617948
#   channels:
#    - bioconda
#   dependencies:
#    - seqprep=1.1
RUN mkdir -p /conda-envs/1f131affa68a486db397c65b7a617948
COPY envs/seqprep.yml /conda-envs/1f131affa68a486db397c65b7a617948/environment.yaml

# Conda environment:
#   source: envs/snpeff.yml
#   prefix: /conda-envs/ae4290c729afa475382c48fb88f731b3
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - snpeff=4.3.1t
RUN mkdir -p /conda-envs/ae4290c729afa475382c48fb88f731b3
COPY envs/snpeff.yml /conda-envs/ae4290c729afa475382c48fb88f731b3/environment.yaml

# Conda environment:
#   source: envs/varscan.yml
#   prefix: /conda-envs/001448ca476ce5fa76876a6c42fc6342
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - varscan=2.4.4
#     - samtools=1.15
#     - tabixpp=1.1.0
RUN mkdir -p /conda-envs/001448ca476ce5fa76876a6c42fc6342
COPY envs/varscan.yml /conda-envs/001448ca476ce5fa76876a6c42fc6342/environment.yaml

# Conda environment:
#   source: envs/vcf2bed.yml
#   prefix: /conda-envs/c5d3e8fb104578aa924e5f6c705c56ab
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bedops=2.4.41
#     - tabixpp=1.1.0
RUN mkdir -p /conda-envs/c5d3e8fb104578aa924e5f6c705c56ab
COPY envs/vcf2bed.yml /conda-envs/c5d3e8fb104578aa924e5f6c705c56ab/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/ecf0f2e646fdd4c5e729de923cf553a6
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - multiqc =1.12
RUN mkdir -p /conda-envs/ecf0f2e646fdd4c5e729de923cf553a6
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/multiqc/environment.yaml /conda-envs/ecf0f2e646fdd4c5e729de923cf553a6/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/3bf9c50269bab91f237f6667f83ab34a --file /conda-envs/3bf9c50269bab91f237f6667f83ab34a/environment.yaml && \
    mamba env create --prefix /conda-envs/e5d860111254c639477a4e8398f307bf --file /conda-envs/e5d860111254c639477a4e8398f307bf/environment.yaml && \
    mamba env create --prefix /conda-envs/f748ab4e3674d4b980a19edf89a39484 --file /conda-envs/f748ab4e3674d4b980a19edf89a39484/environment.yaml && \
    mamba env create --prefix /conda-envs/abfc4906b35726017b76a70cae2c256f --file /conda-envs/abfc4906b35726017b76a70cae2c256f/environment.yaml && \
    mamba env create --prefix /conda-envs/f03630632771bbe472a0c60a72c3b52b --file /conda-envs/f03630632771bbe472a0c60a72c3b52b/environment.yaml && \
    mamba env create --prefix /conda-envs/1f131affa68a486db397c65b7a617948 --file /conda-envs/1f131affa68a486db397c65b7a617948/environment.yaml && \
    mamba env create --prefix /conda-envs/ae4290c729afa475382c48fb88f731b3 --file /conda-envs/ae4290c729afa475382c48fb88f731b3/environment.yaml && \
    mamba env create --prefix /conda-envs/001448ca476ce5fa76876a6c42fc6342 --file /conda-envs/001448ca476ce5fa76876a6c42fc6342/environment.yaml && \
    mamba env create --prefix /conda-envs/c5d3e8fb104578aa924e5f6c705c56ab --file /conda-envs/c5d3e8fb104578aa924e5f6c705c56ab/environment.yaml && \
    mamba env create --prefix /conda-envs/ecf0f2e646fdd4c5e729de923cf553a6 --file /conda-envs/ecf0f2e646fdd4c5e729de923cf553a6/environment.yaml && \
    mamba clean --all -y
