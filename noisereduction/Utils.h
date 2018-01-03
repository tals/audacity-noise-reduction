#pragma once

#include <sndfile.h>
#include <assert.h>

struct SndContext {
    SNDFILE *file;
    SF_INFO info;
};

class SndMmap {
    sf_count_t length;
    int fd;
    char* mmappedData;
    sf_count_t position;
public:
    SF_VIRTUAL_IO interface;
    SndMmap(const char* path);
    SndContext Open();
};

SndContext openAudioFile(const char* path);
