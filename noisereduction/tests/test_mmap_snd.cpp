//
//  test_groundtruth.cpp
//  noisereduction
//
//  Created by Tal Shiri on 1/1/18.
//  Copyright © 2018 org. All rights reserved.
//

#include <stdio.h>
#include <sndfile.h>
#include <cstring>
#include <unistd.h>
#include <cmath>

#include "catch.hpp"
#include "../NoiseReduction.h"
#include "../SndContext.h"
#include "../SndMmap.h"

SndContext openAudioFileRawIO(const char* path) {
    SF_INFO info = { };
    SNDFILE* snd = sf_open(path, SFM_READ, &info);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    return ctx;
}

void checkRead(SndContext& mmapCtx, SndContext& rawCtx, size_t frames, size_t buffSize = 1) {
    for (size_t i = 0; i < frames; i+= buffSize) {
        short bufferMmap[buffSize];
        short bufferRaw[buffSize];

        auto len = sizeof(bufferRaw);

        REQUIRE(sf_readf_short(rawCtx.file, bufferRaw, buffSize));
        REQUIRE(sf_readf_short(mmapCtx.file, bufferMmap, buffSize));
        REQUIRE(memcmp(bufferMmap, bufferRaw, len) == 0);
        REQUIRE(bufferMmap[0] == bufferRaw[0]);
    }
}

TEST_CASE( "MMapped SND", "[MmapSnd]" ) {
    SECTION("Reading everything in one go") {
        SndMmap sndMmaped("dtmf-noise-mono.wav");
        SndContext mmapCtx = sndMmaped.Open();
        SndContext rawCtx = openAudioFileRawIO("dtmf-noise-mono.wav");

        REQUIRE(mmapCtx.file);
        REQUIRE(rawCtx.file);
        REQUIRE(rawCtx.info.frames == mmapCtx.info.frames);
        checkRead(mmapCtx, rawCtx, rawCtx.info.frames);
    }

    SECTION("Seeking") {
        SndMmap sndMmaped("dtmf-noise-mono.wav");
        SndContext mmapCtx = sndMmaped.Open();
        SndContext rawCtx = openAudioFileRawIO("dtmf-noise-mono.wav");

        REQUIRE(mmapCtx.file);
        REQUIRE(rawCtx.file);

        REQUIRE(sf_seek(rawCtx.file, 100, SEEK_SET) == sf_seek(mmapCtx.file, 100, SEEK_SET));
        checkRead(mmapCtx, rawCtx, 10);

        REQUIRE(sf_seek(rawCtx.file, 100, SEEK_CUR) == sf_seek(mmapCtx.file, 100, SEEK_CUR));
        checkRead(mmapCtx, rawCtx, 10);
    }

    SECTION("Read batches") {
        SndMmap sndMmaped("dtmf-noise-mono.wav");
        SndContext mmapCtx = sndMmaped.Open();
        SndContext rawCtx = openAudioFileRawIO("dtmf-noise-mono.wav");

        REQUIRE(mmapCtx.file);
        REQUIRE(rawCtx.file);
        checkRead(mmapCtx, rawCtx, 100, 10);
        checkRead(mmapCtx, rawCtx, 100, 10);
    }
}
