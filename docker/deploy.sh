#!/bin/bash

if [ "$(docker ps -aq -f name=chemical-properties)" != "" ]; then
    echo "removing exited container"
    docker rm -f chemical-properties
fi

docker run -d \
--name chemical-properties \
--restart unless-stopped \
-e ARGS="$*" \
chemical-properties
