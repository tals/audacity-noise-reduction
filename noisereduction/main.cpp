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
#define LOGURU_IMPLEMENTATION 1
#include "loguru.hpp"

int main(int argc, const char * argv[]) {
//    loguru::init(argc, (char**)argv);

    // insert code here...
    SF_INFO info = { 0 };
    SNDFILE* snd = sf_open("/tmp/dtmf.wav", SFM_READ, &info);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    NoiseReduction::Settings settings;
    NoiseReduction reduction(settings, ctx);
    reduction.Process();
    
//    RealFFTf(nullptr, nullptr);
    std::cout << "Hello, World!\n";
    return 0;
}
