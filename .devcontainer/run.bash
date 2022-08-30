#!/bin/bash
. .devcontainer/setup.bash

if docker image inspect $1 > /dev/null 2>&1; then
    :
else
    echo "Image $1 must be created first"
    exit 1
fi

export DISPLAY=:0.0
xhost +local:docker

# --volume "$GITROOT:/app" \
docker run \
    --network=host \
    --env="DISPLAY=$DISPLAY" \
    --env SSH_AUTH_SOCK=$SSH_AUTH_SOCK \
    --volume "$SSH_AUTH_SOCK:$SSH_AUTH_SOCK" \
    --volume "$HOME/.Xauthority:/root/.Xauthority:rw" \
    --volume "/tmp/.X11-unix:/tmp/.X11-unix" \
    --volume "/dev/video0:/dev/video0" \
    --volume "/dev/video1:/dev/video1" \
    --volume "/dev/bus/usb:/dev/bus/usb" \
    -p 8080:80 \
    -it \
    --privileged \
    $1
