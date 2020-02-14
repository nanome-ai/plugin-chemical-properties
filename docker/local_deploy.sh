if [ "$(docker ps -aq -f name=chemical-properties)" != "" ]; then
    # cleanup
    echo "removing exited container"
    docker rm -f chemical-properties
fi

ARGS=$*

docker run -d \
--restart always \
-e ARGS="$ARGS" \
--name chemical-properties chemical-properties
