#!/bin/sh
mkdir -p build
cd build
cmake ..
cmake --build . --target test_runner
:: TODO: clean up: better prepend this path in source code and run directly from build directory
cd ../samples
../build/test_runner
