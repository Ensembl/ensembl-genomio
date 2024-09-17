FROM python:3.10-slim

ARG GENOMIO_VERSION=v1.3.1

RUN apt-get update \
    && apt-get -y upgrade \
    && apt-get -y install \
        default-libmysqlclient-dev \
        pkg-config \
        gcc \
        pip \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir --upgrade pip && pip install --no-cache-dir ensembl-genomio mysqlclient

# Testing dependencies 
CMD ["python", "-c", "'import ensembl.io.genomio'"]
RUN <<EOF
# Test entry point module:
assembly_download --help
EOF

LABEL base.image="ensembl-genomio:${GENOMIO_VERSION}"
LABEL version="$GENOMIO_VERSION"
LABEL software="Ensembl-genomio pipelines"
LABEL software.version="$GENOMIO_VERSION"
LABEL about.summary="A lighweight container to run ensembl-genomio module scripts."
LABEL about.home="https://github.com/Ensembl/ensembl-genomio"
LABEL about.documentation="https://ensembl.github.io/ensembl-genomio"
LABEL license="https://www.apache.org/licenses/LICENSE-2.0.txt"
LABEL mantainer="Ensembl-Metazoa"
LABEL mantainer.email="ensembl-metazoa@ebi.ac.uk"
