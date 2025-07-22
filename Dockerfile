# Interactive environment with scripts and extra dependencies
FROM --platform=linux/amd64 continuumio/miniconda3:24.11.1-0

LABEL org.opencontainers.image.source=https://github.com/formbio/laava

COPY apt-packages.txt /tmp/apt-packages.txt
RUN apt-get update \
    && apt-get install -y findutils \
    && xargs -a /tmp/packages.txt apt-get install -y --no-install-recommends \
    && apt-get remove -y fonts-urw-base35 libgs9 libgs9-common libjbig2dec0 poppler-data \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* /tmp/apt-packages.txt

# Install directly into 'base' conda environment
COPY conda_env.yml ./conda_env.yml
RUN conda env update -v -n base -f conda_env.yml
RUN conda install conda-forge::awscli=2.24.2

# Executable scripts
RUN mkdir -p /opt/laava
RUN chmod +x /opt/laava
COPY src/* /opt/laava/
RUN chmod -R 777 /opt/laava
RUN chmod +x /opt/laava/*.py /opt/laava/*.R /opt/laava/*.sh
ENV PATH="/opt/laava:$PATH"
ENV PYTHONPATH="/opt/laava:$PYTHONPATH"

WORKDIR /data/

CMD ["bash"]
