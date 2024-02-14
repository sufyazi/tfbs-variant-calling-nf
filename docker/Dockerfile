FROM --platform=linux/x86_64 mambaorg/micromamba
USER root
RUN apt-get update && apt-get -y install sudo

LABEL maintainer="Suffian Azizan"

# RUN micromamba config set extract_threads 1

COPY --chown=$MAMBA_USER:$MAMBA_USER env.yaml /tmp/env.yaml

RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
# RUN pip install --no-cache polars pyspark

ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"