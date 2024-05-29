FROM ubuntu:24.04 AS basebuilder

WORKDIR /usr/local/bin/

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install \
        unzip \
        procps \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean


FROM  ubuntu:24.04
COPY --from=basebuilder /usr/local/bin/* /usr/local/bin/
COPY linux-amd64.cli.package.zip /usr/local/bin/
COPY jq /usr/local/bin/

WORKDIR /usr/local/bin/

ARG datasets_version

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


RUN unzip linux-amd64.cli.package.zip \
    && rm linux-amd64.cli.package.zip \
    && chmod +x datasets dataformat \
    && chmod +x jq

CMD ["datasets --version"]
