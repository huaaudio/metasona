// SPDX-License-Identifier: Apache-2.0 AND BSD-3-Clause
// MetaSona adaptation: Jiahua Zhang, September 2026.
/************************************************************************/
/*  Tonality calculation according to Aures (1985)                      */
/*                                                                      */
/*  Ported from Tonality_Aures1985.m (SQAT/PA archive).                 */
/************************************************************************/

#ifndef HEADER_TONALITY_AURES1985
#define HEADER_TONALITY_AURES1985

#ifdef __cplusplus
extern "C"
{
#endif

#if defined(_WIN32) && !defined(TONALITY_AURES1985_STATIC)
#  if defined(TONALITY_AURES1985_BUILD)
#    define TONALITY_AURES1985_API __declspec(dllexport)
#  else
#    define TONALITY_AURES1985_API
#  endif
#else
#  define TONALITY_AURES1985_API
#endif

/* Error codes */
enum _TonalityAures1985ErrorCodes
{
    TonalityAures1985ErrorOutputTooSmall    = -1,
    TonalityAures1985ErrorMemoryAlloc       = -2,
    TonalityAures1985ErrorSignalTooShort    = -3,
    TonalityAures1985ErrorInvalidSampleRate = -4,
    TonalityAures1985ErrorInvalidArgument   = -5,
    TonalityAures1985ErrorFFTPlan           = -6,
    TonalityAures1985ErrorLoudnessFailed    = -7,
    TonalityAures1985ErrorNumerical         = -8
};

/**
 * Compute time-varying Aures tonality from a mono pressure signal.
 *
 * @param signal                Input audio samples (Pa).
 * @param numSamples            Number of samples in signal.
 * @param sampleRate            Sampling rate (Hz). Must be 48 kHz in MetaSona.
 * @param soundField            ISO 532-1 sound field: 0 free, 1 diffuse.
 * @param timeSkip              Seconds to skip for summary statistics.
 * @param outTonality           Output instantaneous tonality (t.u.).
 * @param outTonalWeighting     Output tonal weighting. May be NULL.
 * @param outLoudnessWeighting  Output loudness weighting. May be NULL.
 * @param outTime               Output frame start times in seconds. May be NULL.
 * @param pNumFrames            In: output capacity. Out: frames written.
 * @param outStats              Optional 5 values: mean, std, max, min, p5
 *                              after timeSkip. p5 is the value exceeded
 *                              during 5% of retained frames.
 * @return                      Frames written, or a negative error code.
 */
TONALITY_AURES1985_API int tonality_aures1985(
    const double *signal,
    int           numSamples,
    double        sampleRate,
    int           soundField,
    double        timeSkip,
    double       *outTonality,
    double       *outTonalWeighting,
    double       *outLoudnessWeighting,
    double       *outTime,
    int          *pNumFrames,
    double       *outStats
);

#ifdef __cplusplus
}
#endif

#endif /* HEADER_TONALITY_AURES1985 */
