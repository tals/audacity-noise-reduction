#pragma once

#include "SndContext.h"
#include "InputTrack.h"
#include "OutputTrack.h"

namespace TrackUtils
{
    SndContext openAudioFile(const char* path);
    std::vector<InputTrack> readTracksFromContext(const SndContext& ctx, size_t t0 = 0, size_t t1 = 0);
    InputTrack readOneTrackFromContext(const SndContext &ctx, int channel, size_t t0 = 0, size_t t1 = 0);
    void writeTracksToFile(const char* path, const std::vector<OutputTrack> &tracks, int channels, int sampleRate);
};
