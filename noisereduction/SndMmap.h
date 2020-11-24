#pragma once

#include "SndContext.h"

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
