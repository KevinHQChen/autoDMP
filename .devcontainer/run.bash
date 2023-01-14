#!/bin/bash
. .devcontainer/setup.bash

if docker image inspect $1 > /dev/null 2>&1; then
    :
else
    echo "Image $1 must be created first"
    exit 1
fi

container_name=$(echo $1 | sed -E 's/[^a-zA-Z0-9_.-]//g')

export DISPLAY=:0.0
xhost +local:docker

# check if piezopump exists
if test -e /dev/ttyACM0; then
    export TTYACM0=/dev/ttyACM0
else
    export TTYACM0=/dev/null
fi

# --volume "$GITROOT:/app" \
docker run \
    --name $container_name \
    --network=host \
    --env="DISPLAY=$DISPLAY" \
    --env SSH_AUTH_SOCK=$SSH_AUTH_SOCK \
    --volume "$SSH_AUTH_SOCK:$SSH_AUTH_SOCK" \
    --volume "$HOME/.Xauthority:/root/.Xauthority:rw" \
    --volume "/tmp/.X11-unix:/tmp/.X11-unix" \
    --volume "/dev/video0:/dev/video0" \
    --volume "/dev/video1:/dev/video1" \
    --volume "/dev/bus/usb:/dev/bus/usb" \
    --volume $TTYACM0:$TTYACM0 \
    -p 8080:80 \
    -it \
    --privileged \
    --rm \
    $1
