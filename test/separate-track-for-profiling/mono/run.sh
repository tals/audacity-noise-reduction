#!/bin/sh
input=./square-with-noise.wav
output=./square-with-noise-output.wav
profile=./noise-profile.wav
gain=24
sensitivity=12
smoothing=0
../../../build/noisereduction_driver -i $input -o $output -p $profile --noiseGain $gain --sensitivity $sensitivity --smoothing $smoothing
