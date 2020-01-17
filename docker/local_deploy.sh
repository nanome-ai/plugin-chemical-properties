if [ "$(docker ps -aq -f name=chemical-properties)" != "" ]; then
        # cleanup
        echo "removing exited container"
        docker rm -f chemical-properties
fi

if [ "$1" != "" ]; then
    echo "Using specified plugin server: $1"
    docker run -d \
    -p 8888:8888 \
    -e PLUGIN_SERVER=$1 \
    --name chemical-properties chemical-properties
else
    echo "Using default plugin server: plugins.nanome.ai"
    docker run -d \
    -p 8888:8888 \
    --name chemical-properties chemical-properties
fi