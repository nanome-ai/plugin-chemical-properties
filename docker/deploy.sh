#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

container_name="chemical-properties"
existing=$(docker ps -aq -f name=$container_name)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name $container_name \
--restart unless-stopped \
-h $(hostname)-$container_name \
-e ARGS="$*" \
$container_name
