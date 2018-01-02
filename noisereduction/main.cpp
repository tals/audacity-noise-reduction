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

int main(int argc, const char * argv[]) {
   loguru::g_stderr_verbosity = 9;
   loguru::init(argc, (char**)argv);

    // insert code here...
    SF_INFO info = { 0 };
    SNDFILE* snd = sf_open("/Users/tal/dev/lyrebird/noisereduction/samples/dtmf-noise.wav", SFM_READ, &info);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    NoiseReduction::Settings settings;
    NoiseReduction reduction(settings, ctx);
    reduction.ProfileNoise(3400, 6000);
    reduction.ReduceNoise();
    
    
    return 0;
}
