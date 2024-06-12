FROM ubuntu:24.04 AS basebuilder

ARG datasets_url
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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && echo sh get_datasets_latest.sh ${datasets_url} \
    && sh get_datasets_latest.sh ${datasets_url}

FROM  ubuntu:24.04
COPY --from=basebuilder /usr/local/bin/* /usr/local/bin/
ARG datasets_version="latest"
WORKDIR /usr/local/bin/

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
