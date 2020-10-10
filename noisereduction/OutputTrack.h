#pragma once

#include "Types.h"

class OutputTrack
{
public:
    OutputTrack();
    void Append(float* buffer, size_t length);
    const FloatVector& Buffer() const { return mBuffer; }
    size_t Length() const { return mLength; }
    void SetEnd(size_t newLength);
private:
    FloatVector mBuffer;
    size_t mLength;
};
