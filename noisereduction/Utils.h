#pragma once

#include <sndfile.h>
#include <assert.h>

#define WINDOWS (defined(WIN32) || defined(_WIN32) || defined(__WIN32))

struct SndContext {
    SNDFILE *file;
    SF_INFO info;
};

#ifndef WINDOWS
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
#endif

SndContext openAudioFile(const char* path);
