#include <iostream>
#include "NoiseReduction.h"
#include "InputTrack.h"
#include "OutputTrack.h"
#include "TrackUtils.h"
#include "loguru.hpp"
#include "cxxopts.hpp"

int main(int argc, char * argv[]) {
    cxxopts::Options options("NoiseReduction Test", "Driver for Noise Reduction");
    options.add_options()
        ("i,input", "Input file (required)", cxxopts::value<std::string>())
        ("o,output", "Output file (required)", cxxopts::value<std::string>())
        ("p,profile", "File for noise profiling (must contain only noise)", cxxopts::value<std::string>())
        ("t0", "Start sample for noise profiling", cxxopts::value<size_t>())
        ("t1", "End sample for noise profiling", cxxopts::value<size_t>())
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

    SndContext inputCtx = TrackUtils::openAudioFile(result["input"].as<std::string>().c_str());

    NoiseReduction::Settings settings;
    settings.mNewSensitivity = result["sensitivity"].as<float>();
    settings.mFreqSmoothingBands = result["smoothing"].as<int>();
    settings.mNoiseGain = result["noiseGain"].as<float>();

    NoiseReduction reduction(settings, inputCtx.info.samplerate);

    std::cout << "Profiling noise..." << std::endl;

    std::vector<InputTrack> profileTracks;

    if (result.count("profile")) {
        // use separate track for profiling noise
        SndContext profileCtx = TrackUtils::openAudioFile(result["profile"].as<std::string>().c_str());
        profileTracks = TrackUtils::readTracksFromContext(profileCtx);
    } else {
        // use time indices for profiling noise
        size_t t0 = 0;
        size_t t1 = inputCtx.info.frames;
        if (result["t0"].count()) {
            t0 = result["t0"].as<size_t>();
        }

        if (result["t1"].count()) {
            t1 = result["t1"].as<size_t>();
        }
        
        profileTracks = TrackUtils::readTracksFromContext(inputCtx, t0, t1);
    }

    for (auto& profileTrack : profileTracks) {
        reduction.ProfileNoise(profileTrack);
    }

    std::cout << "Reducing noise..." << std::endl;
    std::vector<InputTrack> inputTracks = TrackUtils::readTracksFromContext(inputCtx);
    std::vector<OutputTrack> outputTracks;
    for (auto& inputTrack : inputTracks) {
        OutputTrack outputTrack;
        reduction.ReduceNoise(inputTrack, outputTrack);
        outputTracks.push_back(outputTrack);
    }
    
    const char* outputPath = result["output"].as<std::string>().c_str();
    TrackUtils::writeTracksToFile(outputPath, outputTracks, inputCtx.info.channels, inputCtx.info.samplerate);

    return 0;
}
