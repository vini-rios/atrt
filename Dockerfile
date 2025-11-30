FROM rocker/tidyverse:4.3.1

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    libgdal-dev \
    python3 \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install conda for Snakemake
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh
ENV PATH="/opt/conda/bin:$PATH"

# Copy environment files
COPY environment.yml /tmp/environment.yml

# Create conda environment
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda env create -f /tmp/environment.yml

# Make scripts executable
COPY scripts/ scripts/
RUN chmod +x scripts/*.sh

# Default command
CMD ["snakemake", "--cores", "all", "--use-conda"]