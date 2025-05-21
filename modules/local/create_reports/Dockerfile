# Stage 1: Get micromamba binary
FROM mambaorg/micromamba:1.5.10-noble AS micromamba

# Stage 2: Main image based on TeX Live
FROM texlive/texlive:latest

# Setup micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=$MAMBA_ROOT_PREFIX/bin:$PATH
ENV PYTHONNOUSERSITE=1

# Copy micromamba binary
COPY --from=micromamba /bin/micromamba /bin/micromamba

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        wget \
        gnupg \
        && rm -rf /var/lib/apt/lists/*

# Install Python & scientific stack with micromamba
RUN micromamba install -y -p $MAMBA_ROOT_PREFIX -c conda-forge \
    python \
    pip \
    pandas \
    numpy \
    plotly \
    pylatex \
    matplotlib \
    scipy \
    && micromamba clean -a -y

# Install pip-only dependencies
RUN pip install --no-cache-dir kaleido==0.1.0

# Set up workspace
WORKDIR /workspace

# Default shell and command
SHELL ["/bin/bash", "-c"]
CMD ["bash"]
