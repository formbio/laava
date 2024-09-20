# Interactive environment with scripts and extra dependencies
FROM --platform=linux/amd64 continuumio/miniconda3:24.7.1-0
LABEL org.opencontainers.image.source https://github.com/formbio/AAV

RUN apt-get update \
    && apt-get install -y \
        apt-transport-https \
        less \
        minimap2 \
        samtools \
        texlive-latex-extra \
        texlive-latex-recommended \
    && apt-get remove -y \
        fonts-urw-base35 \
        libgs9 \
        libgs9-common \
        libjbig2dec0 \
        poppler-data \
    && apt-get autoremove -y
RUN rm -rf /var/lib/apt/lists/*

# Install directly into 'base' conda environment
COPY laava.conda_env.yml ./conda_env.yml
RUN conda env update -v -n base -f conda_env.yml

# Executable scripts
RUN mkdir -p /opt/laava
RUN chmod 777 /opt/laava/
COPY src/* /opt/laava/
RUN chmod +x /opt/laava/*.py /opt/laava/*.R /opt/laava/*.sh
ENV PATH "/opt/laava:$PATH"

WORKDIR /data/

CMD ["bash"]
