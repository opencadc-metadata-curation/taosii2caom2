ARG OPENCADC_PYTHON_VERSION=3.12
FROM opencadc/astropy:${OPENCADC_PYTHON_VERSION}-slim AS builder

RUN apt-get update --no-install-recommends  && apt-get dist-upgrade -y && \
    apt-get install -y build-essential \
                       git \
                       libcurl4-openssl-dev \
                       libhdf5-dev \
                       libnsl-dev \
                       libssl-dev \
                       zlib1g-dev \
    && rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

WORKDIR /usr/src/app

ARG OPENCADC_BRANCH=main
ARG OPENCADC_REPO=opencadc-metadata-curation

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${OPENCADC_BRANCH} && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${OPENCADC_REPO}/taosii2caom2@${OPENCADC_BRANCH}#egg=taosii2caom2

FROM python:${OPENCADC_PYTHON_VERSION}-slim
WORKDIR /usr/src/app
ARG OPENCADC_PYTHON_VERSION

COPY --from=builder /usr/local/lib/python${OPENCADC_PYTHON_VERSION}/site-packages/ /usr/local/lib/python${OPENCADC_PYTHON_VERSION}/site-packages/
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/.config/* /usr/local/.config/

COPY --from=builder /etc/magic /etc/magic
COPY --from=builder /etc/magic.mime /etc/magic.mime
COPY --from=builder /usr/lib/x86_64-linux-gnu/libmagic* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/file/magic.mgc /usr/lib/file/
COPY --from=builder /usr/share/misc/magic /usr/share/misc/magic
COPY --from=builder /usr/share/misc/magic.mgc /usr/share/misc/magic.mgc
COPY --from=builder /usr/share/file/magic.mgc /usr/share/file/magic.mgc

RUN useradd --create-home --shell /bin/bash cadcops
RUN chown -R cadcops:cadcops /usr/src/app

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

