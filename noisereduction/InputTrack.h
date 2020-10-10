#pragma once

#include "Types.h"

class InputTrack
{
public:
    InputTrack(const FloatVector& buffer);
    const FloatVector& Buffer() const { return mBuffer; }
    size_t Length() const { return mBuffer.size(); }
    size_t Read(float* buffer, size_t length);
private:
    FloatVector mBuffer;
    size_t mPosition;
};
