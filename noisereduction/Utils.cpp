#include "Utils.h"
#include <assert.h>
#include <sndfile.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <algorithm>

SndMmap::SndMmap(const char* path) {
    struct stat st;
    stat(path, &st);
    this->length = st.st_size;
    this->fd = open(path, O_RDONLY, 0);
    assert(this->fd);
    this->mmappedData = (char*)mmap(NULL, this->length, PROT_READ, MAP_PRIVATE, this->fd, 0);
    assert(this->mmappedData);
    this->position = 0;

    interface.get_filelen = [](void * user_data) -> sf_count_t {
        SndMmap* self = static_cast<SndMmap*>(user_data);
        return self->length;
    };

    interface.tell = [](void *user_data) -> sf_count_t {
        SndMmap* self = static_cast<SndMmap*>(user_data);
        return self->position;
    };

    interface.seek = [](sf_count_t offset, int whence, void *user_data) -> sf_count_t {
        SndMmap* self = static_cast<SndMmap*>(user_data);
        switch(whence) {
            case SEEK_SET:
                self->position = offset;
                break;
            case SEEK_CUR:
                self->position += offset;
                break;
            case SEEK_END:
                self->position = self->length - offset;
                break;
        }

        self->position = std::max<sf_count_t>(self->position, 0);
        self->position = std::min(self->position, self->length);
        return self->position;
    };

    interface.read = [](void *ptr, sf_count_t count, void *user_data) -> sf_count_t {
        SndMmap* self = static_cast<SndMmap*>(user_data);

        auto len = std::min(count, self->length - self->position);
        memcpy(ptr, self->mmappedData + self->position, len);
        self->position += len;
        return len;
    };
}

SndContext SndMmap::Open() {
    SF_INFO info = { 0 };
    auto snd = sf_open_virtual(&this->interface, SFM_READ, &info, this);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    return ctx;
}

SndContext openAudioFile(const char* path) {
    SF_INFO info = { 0 };
    // SNDFILE* snd = sf_open(path, SFM_READ, &info);
    SndMmap* mmaped = new SndMmap(path);
    auto snd = sf_open_virtual(&mmaped->interface, SFM_READ, &info, mmaped);
    assert(snd);
    SndContext ctx = {
        .file = snd,
        .info = info,
    };

    return ctx;
}


