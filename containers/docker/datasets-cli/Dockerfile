FROM ubuntu:24.04 AS basebuilder

WORKDIR /usr/local/bin/

ARG datasets_version
ARG datasets_download_url
ARG jq_download_url

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install \
        wget \
        unzip \
        procps \
        ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && wget -q $datasets_download_url \
    && unzip linux-amd64.cli.package.zip \
    && rm linux-amd64.cli.package.zip \
    && chmod +x datasets dataformat \
    && wget $jq_download_url -O jq && \
    && chmod +x jq


FROM  ubuntu:24.04
COPY --from=basebuilder /usr/local/bin/data* /usr/local/bin/
COPY --from=basebuilder /usr/local/bin/jq /usr/local/bin/

# WORKDIR /usr/local/bin/

LABEL base.image="ensembl-genomio-datasets:datasets_version"
LABEL version="0.1"
LABEL software="NCBI's datasets and dataformat"
LABEL software.version="$datasets_version"
LABEL about.summary="NCBI Datasets is a new resource that lets you easily gather data from across NCBI databases."
LABEL about.home="https://github.com/ncbi/datasets/releases/tag/v${datasets_version}"
LABEL about.documentation="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
LABEL license="https://spec.openapis.org/oas/v3.0.3"
LABEL mantainer="Ensembl-Metazoa"
LABEL mantainer.email="ensembl-metazoa@ebi.ac.uk"

CMD ["datasets --version"]