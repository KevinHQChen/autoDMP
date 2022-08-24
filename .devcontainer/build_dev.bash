#!/bin/bash
. .devcontainer/setup.bash

docker build \
    --ssh default=${SSH_AUTH_SOCK} \
    -f ./.devcontainer/Dockerfile.dev \
    --tag ubuntu-dev:prebuild .
