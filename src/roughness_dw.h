// SPDX-License-Identifier: Apache-2.0 AND MIT
// MetaSona adaptation: Jiahua Zhang and Codex, September 2026.
/************************************************************************/
/*  Roughness calculation according to Daniel & Weber (1997)            */
/*  "Psychoacoustical roughness: implementation of an optimized model"  */
/*  Acustica(83), 113-123.                                              */
/*                                                                      */
/*  Reference: 60 dB 1 kHz tone 100% modulated at 70 Hz = 1 asper      */
/************************************************************************/

#ifndef HEADER_ROUGHNESS_DW
#define HEADER_ROUGHNESS_DW

#ifdef __cplusplus
extern "C"
{
#endif

/* Number of critical band channels (0.5 Bark resolution, 47 channels) */
#define ROUGHNESS_N_CHANNELS    47

/* Error codes */
enum _RoughnessErrorCodes
{
    RoughnessErrorOutputTooSmall  = -1,
    RoughnessErrorMemoryAlloc     = -2,
    RoughnessErrorSignalTooShort  = -3,
    RoughnessErrorInvalidSampleRate = -4,
    RoughnessErrorNumerical = -5
};

/**
 * Compute time-varying roughness from a mono audio signal.
 *
 * @param pSignal       Pointer to input audio samples (Pa).
 * @param numSamples    Number of samples in pSignal.
 * @param sampleRate    Sampling rate (Hz). Supported: 44100, 40960, 48000
 *                      (other rates are resampled to 48000 internally).
 * @param outRoughness  Output array for instantaneous roughness (asper).
 *                      Must have space for at least *pNumFrames doubles.
 * @param outSpecRoughness  Output array for time-averaged specific roughness
 *                          (asper/Bark), 47 values (one per 0.5 Bark band).
 *                          May be NULL if not needed.
 * @param pNumFrames    In: maximum number of frames that fit in outRoughness.
 *                      Out: actual number of frames written.
 * @return              Number of frames written, or negative error code.
 */
int roughness_dw(
    const double *pSignal,
    int           numSamples,
    double        sampleRate,
    double       *outRoughness,
    double       *outSpecRoughness,
    int          *pNumFrames
);

#ifdef __cplusplus
}
#endif

#endif /* HEADER_ROUGHNESS_DW */
