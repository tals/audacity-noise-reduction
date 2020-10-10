#include "OutputTrack.h"
#include <assert.h>

OutputTrack::OutputTrack() :
    mLength(0)
{ }

void OutputTrack::Append(float* buffer, size_t length)
{
    mBuffer.insert(mBuffer.end(), buffer, &buffer[length]);
    mLength += length;
}

void OutputTrack::SetEnd(size_t newLength)
{
    assert(newLength <= mLength);
    mLength = newLength;
}
