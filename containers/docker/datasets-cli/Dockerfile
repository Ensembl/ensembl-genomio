FROM ubuntu:24.04 AS basebuilder

WORKDIR /usr/local/bin/

ADD get_datasets_latest.sh /usr/local/bin/

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install \
        wget \
        unzip \
        procps \
        ca-certificates \
        grep=3.1* \
        curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && sh get_datasets_latest.sh

FROM  ubuntu:24.04
COPY --from=basebuilder /usr/local/bin/* /usr/local/bin/

WORKDIR /usr/local/bin/

RUN . /usr/local/bin/datasets_version.sh

ENV datasets_version=$DATASETS_RELEASE

CMD ["datasets --version"]

LABEL base.image="ensembl-genomio-datasets:${datasets_version}"
LABEL version="0.1"
LABEL software="NCBI's datasets and dataformat"
LABEL software.version="$datasets_version"
LABEL about.summary="NCBI Datasets is a new resource that lets you easily gather data from across NCBI databases."
LABEL about.home="https://github.com/ncbi/datasets/releases/tag/v${datasets_version}"
LABEL about.documentation="https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
LABEL license="https://spec.openapis.org/oas/v3.0.3"
LABEL mantainer="Ensembl-Metazoa"
LABEL mantainer.email="ensembl-metazoa@ebi.ac.uk"