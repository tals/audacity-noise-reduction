#!/bin/sh
mkdir -p build
cd build
cmake ..
cmake --build . --target noisereduction_driver
