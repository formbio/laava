FROM ubuntu:22.04

ENV QUALIMAP_VER="2.3"
ENV FASTQC_VER="0.12.1"
#Build Qualimap Environment variable

RUN apt-get update; apt-get install -y autoconf automake build-essential gcc checkinstall; apt-get upgrade; \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --allow-unauthenticated libxml2-dev libssl-dev libcurl4-gnutls-dev wget unzip git default-jre default-jdk r-base pigz parallel bzip2 curl gcc libc6-dev  make zlib1g samtools bedtools fastqc && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN apt-get update; apt-get remove -y ghostscript poppler-data fonts-urw-base35 libjbig2dec0 libpulse0

RUN apt update -qq; DEBIAN_FRONTEND="noninteractive" apt install -y --allow-unauthenticated software-properties-common dirmngr

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu/focal-cran40/"; apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9; add-apt-repository ppa:openjdk-r/ppa

RUN apt update; DEBIAN_FRONTEND="noninteractive" apt install -y --allow-unauthenticated r-base r-base-core r-recommended r-base-dev  && rm -rf /var/lib/apt/lists/*

RUN mkdir -p /opt;
#Install Qualimap
RUN cd /opt && wget https://bitbucket.org/kokonech/qualimap/downloads/qualimap_v${QUALIMAP_VER}.zip
RUN cd /opt && unzip qualimap_v${QUALIMAP_VER}.zip && rm qualimap_v${QUALIMAP_VER}.zip

RUN R -e 'install.packages(c("optparse","XML","BiocManager","remotes","RColorBrewer","tidyverse"), repos = "http://cran.r-project.org")'
RUN Rscript -e 'BiocManager::install(c("NOISeq","Repitools","Rsamtools","rtracklayer","GenomicFeatures"))'

ENV PATH "$PATH:/opt/qualimap_v${QUALIMAP_VER}"
ENV PATH "$PATH:/usr/local/bin"

WORKDIR /data/
