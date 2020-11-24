#pragma once

#include <sndfile.h>

struct SndContext {
    SNDFILE* file;
    SF_INFO info;
};
