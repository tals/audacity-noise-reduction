#include "Utils.h"
#include <assert.h>

SndContext openAudioFile(const char* path) {
    SF_INFO info = { 0 };
    SNDFILE* snd = sf_open(path, SFM_READ, &info);
    assert(snd);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    return ctx;
}


