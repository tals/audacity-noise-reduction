//
//  main.cpp
//  noisereduction
//
//  Created by Tal Shiri on 12/29/17.
//  Copyright © 2017 org. All rights reserved.
//

#include <iostream>
#include "NoiseReduction.h"
#include <sndfile.h>
#include "loguru.hpp"
#include "Utils.h"
int main(int argc, const char * argv[]) {
   loguru::g_stderr_verbosity = 9;
   loguru::init(argc, (char**)argv);

    // insert code here...
    SF_INFO info = { 0 };
    SNDFILE* snd = sf_open("/Users/tal/dev/lyrebird/noisereduction/samples/dtmf-noise-stereo.wav", SFM_READ, &info);
    assert(snd);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };
    NoiseReduction::Settings settings;

    // test settings
    settings.mNewSensitivity = 16;
    settings.mFreqSmoothingBands = 0;
    settings.mNoiseGain = 39;

    NoiseReduction reduction(settings, ctx);
    reduction.ProfileNoise(3400, 6000);
    reduction.ReduceNoise("/tmp/foo.wav");


    return 0;
}
