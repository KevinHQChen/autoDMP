.DEFAULT_GOAL := help   # set default target if no arguments are given to make

.PHONY: help configure build test test_release debug release docs format clean

help:	## Show help.
	@grep -hE '^[A-Za-z0-9_ \-]*?:.*##.*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

configure:	## Configure the build.
	make debug_config

build:	## Build the project.
	make debug

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

#   -DGIT_SHA:STRING=${{ github.sha }} # https://medium.com/@mtiller/using-git-hashes-in-makefile-rules-387a099b9cb
debug_config:
	cmake -S . -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=Debug -DENABLE_DEVELOPER_MODE:BOOL=OFF -DOPT_ENABLE_COVERAGE:BOOL=ON

debug:
	cmake --build ./build --config Debug

release:
	cmake -S ./ -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=Release -DFEATURE_TESTS:BOOL=OFF
	cmake --build ./build --config Release

test:
	cmake -S ./ -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=Debug -DFEATURE_TESTS:BOOL=ON
	cmake --build ./build --config Debug

	(cd build/my_exe/test && ctest -C Debug --output-on-failure)
	(cd build/my_header_lib/test && ctest -C Debug --output-on-failure)
	(cd build/my_lib/test && ctest -C Debug --output-on-failure)

test_release_debug:
	cmake -S ./ -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo -DFEATURE_TESTS:BOOL=ON
	cmake --build ./build --config RelWithDebInfo

	(cd build/my_exe/test && ctest -C RelWithDebInfo --output-on-failure)
	(cd build/my_header_lib/test && ctest -C RelWithDebInfo --output-on-failure)
	(cd build/my_lib/test && ctest -C RelWithDebInfo --output-on-failure)

test_release:
	cmake -S ./ -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=Release -DFEATURE_TESTS:BOOL=ON
	cmake --build ./build --config Release

	(cd build/my_exe/test && ctest -C Release --output-on-failure)
	(cd build/my_header_lib/test && ctest -C Release --output-on-failure)
	(cd build/my_lib/test && ctest -C Release --output-on-failure)

test_install:
	cmake --install ./build --prefix ./build/test_install

coverage:
	make test
	gcovr -j 1 --delete --root ./ --print-summary --xml-pretty --xml coverage.xml ./build --gcov-executable gcov

docs:
	cmake -S ./ -B ./build -G "Ninja Multi-Config" -DCMAKE_BUILD_TYPE:STRING=Debug -DFEATURE_DOCS:BOOL=ON -DFEATURE_TESTS:BOOL=OFF
	cmake --build ./build --target doxygen-docs --config Debug

format: ## Format source code using clang-format.
	git ls-files --exclude-standard | grep -E '\.(cpp|hpp|c|cc|cxx|hxx|ixx)$$' | xargs clang-format -i -style=file

clean:
	rm -rf ./build
