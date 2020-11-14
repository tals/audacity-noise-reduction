#!/bin/sh
input=./square-with-noise.wav
output=./square-with-noise-output.wav
t0=0
t1=44100
gain=24
sensitivity=12
smoothing=0
../../../build/noisereduction_driver -i $input -o $output --t0 $t0 --t1 $t1 --noiseGain $gain --sensitivity $sensitivity --smoothing $smoothing
