# autoDMP

An automatic droplet microfluidics platform.

## Usage

- download and install docker, add current user to docker group so we don't need sudo.
- ensure buildx is used (```docker buildx install```).
- run ```make docker-build``` to build Docker image with relevant build tools.
- run ```make docker-run``` to run Docker container from current image.
- TODO: add further documentation
