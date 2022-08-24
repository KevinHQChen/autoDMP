#!/bin/bash
GITROOT=$(git rev-parse --show-toplevel)
GITSHA=$(git rev-parse HEAD | head -c7)
export GIT_EMAIL=$(git config --get user.email)
export GIT_NAME=$(git config --get user.name)

export DOCKER_BUILDKIT=1
export BUILDKIT_PROGRESS=plain
