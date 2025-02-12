#!/bin/bash

COLLECTION="taosii"
IMAGE="opencadc/${COLLECTION}2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE}

echo "Run image ${IMAGE}"
docker run --rm --name ${COLLECTION}_todo --user $(id -u):$(id -g) -e HOME=/usr/src/app --mount type=bind,src=${PWD},dst=/usr/src/app/ --mount type=bind,src=${PWD}/test_files,dst=/data ${IMAGE} ${COLLECTION}_run || exit $?

date
exit 0
