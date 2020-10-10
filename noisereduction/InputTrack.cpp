#include "InputTrack.h"

InputTrack::InputTrack(const FloatVector& buffer) :
    mBuffer(buffer),
    mPosition(0)
{ }

size_t InputTrack::Read(float* buffer, size_t length)
{
    float *writePtr = buffer;
    size_t totalRead = 0;

    for (size_t i = 0; i < length; i++) {
        if (mPosition >= mBuffer.size()) {
            // reached end of this buffer
            break;
        }

        *writePtr = mBuffer[mPosition];
        mPosition++;
        writePtr++;
        totalRead++;
    }

    return totalRead;
}
