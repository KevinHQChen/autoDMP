#!/bin/bash
. .devcontainer/setup.bash

docker build \
    --ssh default=${SSH_AUTH_SOCK} \
    --build-arg VARIANT=$1 \
    --build-arg GIT_EMAIL \
    --build-arg GIT_NAME \
    -f ./.devcontainer/Dockerfile.ubuntu-cpp \
    --tag ubuntu-cpp:prebuild .
