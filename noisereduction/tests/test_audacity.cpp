//
//  test_groundtruth.cpp
//  noisereduction
//
//  Created by Tal Shiri on 1/1/18.
//  Copyright © 2018 org. All rights reserved.
//

#include <stdio.h>
#include <sndfile.h>
#include <string>
#include <cmath>

#include "catch.hpp"
#include "../TrackUtils.h"
#include "../Utils.h"
#include "../NoiseReduction.h"

std::vector<float> readContext(SndContext& ctx) {
    auto size = ctx.info.frames * ctx.info.channels;
    std::vector<float> buff(size);
    size_t read = sf_read_float(ctx.file, &buff[0], size);
    assert(read == size);

    return buff;
}

void compare(const char* inputPath, const char* groundTruthPath) {
    SndContext inputCtx = openAudioFile(inputPath);
    SndContext groundTruthCtx = openAudioFile(groundTruthPath);

    NoiseReduction::Settings settings;
    settings.mNewSensitivity = 16;
    settings.mFreqSmoothingBands = 0;
    settings.mNoiseGain = 39;

    auto const t0 = 3400;
    auto const t1 = 6000;
    std::vector<InputTrack> profileTracks = TrackUtils::readTracksFromContext(inputCtx, t0, t1);

    NoiseReduction reduction(settings, inputCtx.info.samplerate);
    for (auto& profileTrack : profileTracks) {
        reduction.ProfileNoise(profileTrack);
    }

    std::vector<InputTrack> inputTracks = TrackUtils::readTracksFromContext(inputCtx);
    std::vector<OutputTrack> outputTracks;
    for (auto& inputTrack : inputTracks) {
        OutputTrack outputTrack;
        reduction.ReduceNoise(inputTrack, outputTrack);
        outputTracks.push_back(outputTrack);
    }

    TrackUtils::writeTracksToFile("/tmp/processed.wav", outputTracks, inputCtx.info.channels, inputCtx.info.samplerate);

    SndContext processedCtx = openAudioFile("/tmp/processed.wav");
    REQUIRE(processedCtx.info.frames == groundTruthCtx.info.frames);
    REQUIRE(processedCtx.info.channels == groundTruthCtx.info.channels);

    auto groundTruthData = readContext(groundTruthCtx);
    auto processedData = readContext(processedCtx);

    const float EPS = 0.0005;
    for (int i = 0; i < processedData.size(); i++) {
        auto gtSample = groundTruthData[i];
        auto processedSample = processedData[i];
        auto delta = std::fabs(processedSample - gtSample);
        REQUIRE(delta < EPS);
    }
}

TEST_CASE( "Noise Reduction", "[NoiseReduction]" ) {
    SECTION( "mono track" ) {
        compare("dtmf-noise-mono.wav", "dtmf-noise-mono-audacity-gain-39-sensitivity-16-smooth-0.wav");
    }

    SECTION( "stereo track" ) {
        compare("dtmf-noise-stereo.wav", "dtmf-noise-stereo-audacity-gain-39-sensitivity-16-smooth-0.wav");
    }
}
