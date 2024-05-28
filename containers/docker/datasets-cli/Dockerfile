FROM ubuntu:24.04 AS basebuilder

WORKDIR /usr/local/bin/

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install \
        wget \
        unzip \
        procps \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean


FROM  ubuntu:24.04
COPY --from=basebuilder /usr/local/bin/* /usr/local/bin/
COPY get_datasets_latest.sh /usr/local/bin/

WORKDIR /usr/local/bin/
RUN source /usr/local/bin/get_datasets_latest.sh

ENV DATASETS_VERSION=$DATASETS_RELEASE
ENV DATASETS_DOWNLOAD_URL=$DATASETS_LATEST
ENV JQ_DOWNLOAD_URL=$JQ_LATEST

LABEL base.image="ensembl-genomio-datasets:latest"
LABEL version="0.1"
LABEL software="NCBI's datasets and dataformat"
LABEL software.version="${DATASETS_VERSION}"
LABEL about.summary="NCBI Datasets is a new resource that lets you easily gather data from across NCBI databases."
LABEL about.home="https://github.com/ncbi/datasets/releases/tag/v${DATASETS_VERSION}"
LABEL about.documentation="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
LABEL license="https://spec.openapis.org/oas/v3.0.3"
LABEL mantainer="Ensembl-Metazoa"
LABEL mantainer.email="ensembl-metazoa@ebi.ac.uk"

RUN wget -q ${DATASETS_DOWNLOAD_URL} \
    && linux-amd64.cli.package.zip \
    && rm linux-amd64.cli.package.zip \
    && chmod +x datasets dataformat \
    && wget ${JQ_DOWNLOAD_URL} -O jq && \
    && chmod +x jq

CMD ["datasets --version"]