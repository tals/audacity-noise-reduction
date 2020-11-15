#!/bin/sh
mkdir -p build
cd build
cmake ..
cmake --build . --target test_runner
cd ../samples
../build/test_runner
