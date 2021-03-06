#include "TrackUtils.h"

#include <memory>
#include <stdexcept>

#include "MemoryX.h"
#include "loguru.hpp"

#ifdef USE_MMAP
#include "SndMmap.h"
#endif

SndContext TrackUtils::openAudioFile(const char* path) {
    SF_INFO info = {};
#ifdef USE_MMAP
    SndMmap* mmaped = new SndMmap(path);
    SNDFILE* snd = sf_open_virtual(&mmaped->interface, SFM_READ, &info, mmaped);
#else
    // no mmap support, so fall back to simple solution
    SNDFILE* snd = sf_open(path, SFM_READ, &info);
#endif
    assert(snd);
    SndContext ctx = {};
    ctx.file = snd;
    ctx.info = info;

    return ctx;
}

std::vector<InputTrack> TrackUtils::readTracksFromContext(const SndContext& ctx, size_t t0/* = 0*/, size_t t1/* = 0*/)
{
    std::vector<InputTrack> tracks;

    for (int channel = 0; channel < ctx.info.channels; channel++) {
        LOG_F(INFO, "Reading channel %d", channel);
        InputTrack track = readOneTrackFromContext(ctx, channel, t0, t1);
        tracks.push_back(track);
    }

    return tracks;
}

InputTrack TrackUtils::readOneTrackFromContext(const SndContext &ctx, int channel, size_t t0/* = 0*/, size_t t1/* = 0*/)
{
    // if t1 is undefined, read full track
    if (t1 == 0)
        t1 = (size_t)ctx.info.frames;

    const size_t frameCount = t1 - t0;
    FloatVector buffer(frameCount);

    sf_seek(ctx.file, t0, SEEK_SET);
    float* writePtr = &buffer[0];

    auto frameBuffer = std::make_unique<float[]>(ctx.info.channels);

    // can only read full frames from libsnd.
    // will probably be a lot faster to read the whole thing and then split it
    for (size_t frame = t0; frame < t1; frame++) {
        size_t read = sf_readf_float(ctx.file, frameBuffer.get(), 1);

        if (read == 0) {
            LOG_F(WARNING, "Couldn't read data from input file");
            break;
        }

        *writePtr = frameBuffer[channel];
        writePtr++;
    }

    return InputTrack(buffer);
}

void TrackUtils::writeTracksToFile(const char* path, const std::vector<OutputTrack> &tracks, int channels, int sampleRate)
{
    if (tracks.empty()) {
        LOG_F(WARNING, "Tracks are empty");
        return;
    }

    SF_INFO info = { };
    info.channels = channels;
    info.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
    info.samplerate = sampleRate;

    SNDFILE* sf = sf_open(path, SFM_WRITE, &info);
    if (!sf) {
        throw std::runtime_error("Cannot open output file");
    }

    size_t frameCount = tracks[0].Length();
    LOG_F(INFO, "Writing %zd frames to file", frameCount);

    // copy audio to buffer so that channels are interleaved
    const size_t bufferFrames = 1024;
    const size_t bufferSize = channels * bufferFrames;
    auto buffer = std::make_unique<float[]>(bufferSize);

    int frames = 0;

    for (size_t frame = 0; frame < frameCount; frame++) {
        for (int channel = 0; channel < channels; channel++) {
            buffer[frames * channels + channel] = tracks[channel].Buffer()[frame];
        }

        frames++;

        if (frames == bufferFrames) {
            sf_count_t written = sf_writef_float(sf, buffer.get(), bufferFrames);
            assert(written > 0);
            frames = 0;
        }
    }

    // write remaining frames left in buffer
    if (frames > 0) {
        sf_count_t written = sf_writef_float(sf, buffer.get(), frames);
        assert(written > 0);
    }

    sf_close(sf);
}
