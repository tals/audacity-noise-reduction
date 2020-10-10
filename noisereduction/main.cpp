#include <iostream>
#include "NoiseReduction.h"
#include "InputTrack.h"
#include "OutputTrack.h"
#include "TrackUtils.h"
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

    SndContext ctx = openAudioFile(result["input"].as<std::string>().c_str());

    NoiseReduction::Settings settings;
    settings.mNewSensitivity = result["sensitivity"].as<float>();
    settings.mFreqSmoothingBands = result["smoothing"].as<int>();
    settings.mNoiseGain = result["noiseGain"].as<float>();

    NoiseReduction reduction(settings, ctx.info.samplerate);

    size_t t0 = 0;
    size_t t1 = ctx.info.frames;
    if (result["t0"].count()) {
        t0 = result["t0"].as<size_t>();
    }

    if (result["t1"].count()) {
        t1 = result["t1"].as<size_t>();
    }

    std::cout << "Profiling noise..." << std::endl;
    std::vector<InputTrack> profileTracks = TrackUtils::readTracksFromContext(ctx, t0, t1);
    for (auto& profileTrack : profileTracks) {
        reduction.ProfileNoise(profileTrack);
    }

    std::cout << "Denoising..." << std::endl;
    std::vector<InputTrack> inputTracks = TrackUtils::readTracksFromContext(ctx);
    std::vector<OutputTrack> outputTracks;
    for (auto& inputTrack : inputTracks) {
        OutputTrack outputTrack;
        reduction.ReduceNoise(inputTrack, outputTrack);
        outputTracks.push_back(outputTrack);
    }
    
    const char* outputPath = result["output"].as<std::string>().c_str();
    TrackUtils::writeTracksToFile(outputPath, outputTracks, ctx.info.channels, ctx.info.samplerate);

    return 0;
}
