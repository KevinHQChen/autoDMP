CMAKE_BUILD_TYPE ?= Release	# set the default build type to Release if none was specified
CMAKE_PRESET ?= Multi-Config		# use the Multi-Config preset by default (configures all build types at once)
BUILD_TESTS ?= OFF		# disable tests by default

.DEFAULT_GOAL := help   	# set default target if no arguments are given to make

# these targets will always be executed when called, even if a file with the same name exists
.PHONY: help configure build debug release relwithdebinfo test test_debug test_release test_relwithdebinfo test_install coverage docs build_docs format clean docker-build docker-build-dev docker-build-deps docker-run-prebuild docker-run

help:					## Show help.
	@grep -hE '^[A-Za-z0-9_ \-]*?:.*##.*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

configure:				## Configure the build (default build type: Release).
	cmake --preset=$(CMAKE_PRESET)

build: configure		## Build the project (default build type: Release).
	cmake --build ./build --config $(CMAKE_BUILD_TYPE)

#   -DGIT_SHA:STRING=${{ github.sha }} # https://medium.com/@mtiller/using-git-hashes-in-makefile-rules-387a099b9cb
debug: 					## Build the project in debug mode.
	$(MAKE) CMAKE_BUILD_TYPE=Debug build

release: 				## Build the project in release mode.
	$(MAKE) CMAKE_BUILD_TYPE=Release build

relwithdebinfo: 		## Build the project in release with debug info mode.
	$(MAKE) CMAKE_BUILD_TYPE=RelWithDebInfo build

test: 					## Run tests (default build type: Release).
	$(MAKE) BUILD_TESTS=ON build
	(cd build/autoDMP/test && ctest -C $(CMAKE_BUILD_TYPE) --output-on-failure)
	(cd build/util/test && ctest -C $(CMAKE_BUILD_TYPE) --output-on-failure)
	(cd build/gui/test && ctest -C $(CMAKE_BUILD_TYPE) --output-on-failure)
	(cd build/cam/test && ctest -C $(CMAKE_BUILD_TYPE)  --output-on-failure)
	(cd build/ctrl/test && ctest -C $(CMAKE_BUILD_TYPE) --output-on-failure)

test_debug: 			## Run tests in debug mode.
	$(MAKE) CMAKE_BUILD_TYPE=Debug BUILD_TESTS=ON test

test_release: 			## Run tests in release mode.
	$(MAKE) CMAKE_BUILD_TYPE=Release BUILD_TESTS=ON test

test_relwithdebinfo: 	## Run tests in release with debug info mode.
	$(MAKE) CMAKE_BUILD_TYPE=RelWithDebInfo BUILD_TESTS=ON test

test_install:
	cmake --install ./build --prefix ./build/test_install

coverage: 				## Run tests and generate coverage report.
	make test
	gcovr -j 1 --delete --root ./ --print-summary --xml-pretty --xml coverage.xml ./build --gcov-executable gcov

docs: 					## Build the documentation.
	$(MAKE) CMAKE_BUILD_TYPE=Debug CMAKE_PRESET=docs configure
	$(MAKE) CMAKE_BUILD_TYPE=Debug build_docs

build_docs:
	cmake --build ./build --config $(CMAKE_BUILD_TYPE) --target docs

format:
	git ls-files --exclude-standard | grep -E '\.(cpp|hpp|c|cc|cxx|hxx|ixx)$$' | xargs clang-format -i -style=file

clean:
	rm -rf ./build
	rm -rf ./ctrl/scripts/venv

docker-build:	## Build ubuntu-cpp:prebuild image (build tools).
	./.devcontainer/build_deps.bash "ubuntu:20.04"

docker-build-dev:	## Build ubuntu-dev:prebuild image (custom dev tools).
	./.devcontainer/build_dev.bash

docker-build-deps:	## Build ubuntu-cpp:prebuild image (both dev and build tools).
	./.devcontainer/build_deps.bash "ubuntu-dev:prebuild"

docker-run-prebuild:	## Run ubuntu-cpp:prebuild container from current image.
	./.devcontainer/run.bash "ubuntu-cpp:prebuild"

docker-run:	## Run ubuntu-cpp:latest container from current image.
	./.devcontainer/run.bash "ubuntu-cpp:latest"

docker-commit:	## Commit current container as ubuntu-cpp:latest.
	docker commit ubuntu-cpplatest ubuntu-cpp:latest
