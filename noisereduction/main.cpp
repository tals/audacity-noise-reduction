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
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    cxxopts::Options options("NoiseReduction Test", "Driver for Noise Reduction");
    options.add_options()
        ("i,input", "Input file (required)", cxxopts::value<std::string>())
        ("o,output", "Output file (required)", cxxopts::value<std::string>())
        ("t0", "Start sample for noise profiling", cxxopts::value<size_t>())
        ("t1", "end sample for noise profiling", cxxopts::value<size_t>())
        ("noiseGain", "Noise Gain (dB)", cxxopts::value<float>()->default_value("12"))
        ("sensitivity", "Sensitivity", cxxopts::value<float>()->default_value("6.0"))
        ("smoothing", "Frequency Smoothing (bands)", cxxopts::value<int>()->default_value("0"))
        ("verbose", "Verbose Output");

    options.parse_positional(std::vector<std::string>{"input", "output"});

    int unparsedArgc = argc;
    auto result = options.parse(unparsedArgc, argv);
    if (!result.count("input") || !result.count("output")) {
        auto help = options.help();
        std::cout << help << std::endl;
        return 1;
    }

    loguru::g_stderr_verbosity = loguru::Verbosity_ERROR;
    if (result["verbose"].count()) {
        loguru::g_stderr_verbosity = loguru::Verbosity_1;
    }

    loguru::init(argc, argv);


    std::cout << "Processing " << result["input"].as<std::string>() << " -> " << result["output"].as<std::string>() << std::endl;

    auto ctx = openAudioFile(result["input"].as<std::string>().c_str());

    NoiseReduction::Settings settings;
    settings.mNewSensitivity = result["sensitivity"].as<float>();
    settings.mFreqSmoothingBands = result["smoothing"].as<int>();
    settings.mNoiseGain = result["noiseGain"].as<float>();
    NoiseReduction reduction(settings, ctx);
    auto t0 = 0;
    auto t1 = ctx.info.frames;

    if (result["t0"].count()) {
        t0 = result["t0"].as<size_t>();
    }

    if (result["t1"].count()) {
        t1 = result["t1"].as<size_t>();
    }

    std::cout << "Profiling noise..." << std::endl;
    reduction.ProfileNoise(t0, t1);
    std::cout << "Denoising..." << std::endl;
    reduction.ReduceNoise(result["output"].as<std::string>().c_str());


    return 0;
}
