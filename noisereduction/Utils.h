#pragma once

#include <sndfile.h>
#include <assert.h>

struct SndContext {
    SNDFILE *file;
    SF_INFO info;
};

SndContext openAudioFile(const char* path);
