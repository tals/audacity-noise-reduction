#pragma once
#include <sndfile.h>
#include <memory>
#include "Utils.h"
#define DB_TO_LINEAR(x) (pow(10.0, (x) / 20.0))
#define LINEAR_TO_DB(x) (20.0 * log10(x))

typedef char *samplePtr;

class NoiseReductionWorker;
class Statistics;
class NoiseReduction {
public:
    struct Settings {
        Settings();

        size_t WindowSize() const { return 1u << (3 + mWindowSizeChoice); }
        unsigned StepsPerWindow() const { return 1u << (1 + mStepsPerWindowChoice); }
        bool       mDoProfile;
        double     mNewSensitivity;   // - log10 of a probability... yeah.
        double     mFreqSmoothingBands; // really an integer
        double     mNoiseGain;         // in dB, positive
        double     mAttackTime;        // in secs
        double     mReleaseTime;       // in secs

        // Advanced:
        double     mOldSensitivity;    // in dB, plus or minus

        // Basic:
        int        mNoiseReductionChoice;

        // Advanced:
        int        mWindowTypes;
        int        mWindowSizeChoice;
        int        mStepsPerWindowChoice;
        int        mMethod;
    };

    NoiseReduction(NoiseReduction::Settings& settings, SndContext& ctx);
    ~NoiseReduction() {};
    void ProfileNoise(size_t t0, size_t t1);
    void ReduceNoise(const char* outputPath, size_t t0, size_t t1);
    void ReduceNoise(const char* outputPath);
private:
    std::auto_ptr<Statistics> mStatistics;
    NoiseReduction::Settings mSettings;
    SndContext& mCtx;

};
