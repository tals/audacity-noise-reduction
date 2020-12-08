#!/bin/sh
mkdir -p build
cd build
cmake ..
cmake --build . --target test_runner
mkdir -p temp
./test_runner
