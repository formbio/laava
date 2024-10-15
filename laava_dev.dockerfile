# Development environment for running the scripts, no scripts, extra dependencies
FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0

# Set the container's timezone to match this local machine
RUN ln -snf /usr/share/zoneinfo/$CONTAINER_TIMEZONE /etc/localtime && echo $CONTAINER_TIMEZONE > /etc/timezone
# Silence a debconf warning
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        build-essential \
        git \
        less \
        minimap2 \
        samtools \
        texlive-latex-extra \
        texlive-latex-recommended \
        vim \
    && apt-get remove -y \
        fonts-urw-base35 \
        libgs9 \
        libgs9-common \
        libjbig2dec0 \
        poppler-data \
    && apt-get autoremove -y

# Install directly into 'base' conda environment
COPY laava_dev.conda_env.yml ./conda_env.yml
RUN conda env update -v -n base -f conda_env.yml

WORKDIR /data/

CMD ["bash"]
