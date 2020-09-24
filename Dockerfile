FROM opencadc/astropy:3.8-slim

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    git
    
RUN pip install cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    caom2utils \
    ftputil \
    importlib-metadata \
    pytz \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

RUN apt-get install -y libhdf5-dev

RUN pip install h5py

ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2pipe.git && \
  pip install ./caom2pipe
  
RUN git clone https://github.com/${OPENCADC_REPO}/taosii2caom2.git && \
  cp ./taosii2caom2/scripts/config.yml / && \
  cp ./taosii2caom2/scripts/docker-entrypoint.sh / && \
  pip install ./taosii2caom2

ENTRYPOINT ["/docker-entrypoint.sh"]
