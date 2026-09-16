// SPDX-License-Identifier: Apache-2.0
// MetaSona authors: Jiahua Zhang and Codex, September 2026.

#ifndef METASONA_METASONA_H
#define METASONA_METASONA_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32) || defined(__CYGWIN__)
#  define MS_CALL __cdecl
#  if defined(MS_BUILDING_LIBRARY)
#    define MS_API __declspec(dllexport)
#  elif defined(MS_STATIC)
#    define MS_API
#  else
#    define MS_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define MS_CALL
#  define MS_API __attribute__((visibility("default")))
#else
#  define MS_CALL
#  define MS_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define MS_ABI_VERSION 1u
#define MS_SAMPLE_RATE_HZ 48000u
#define MS_THIRD_OCTAVE_BANDS 28u
#define MS_BARK_BANDS 240u

typedef int32_t ms_status;
enum ms_status_value {
    MS_OK = 0,
    MS_ERROR_NULL_POINTER = 1,
    MS_ERROR_INVALID_ARGUMENT = 2,
    MS_ERROR_INVALID_SAMPLE_RATE = 3,
    MS_ERROR_INVALID_SOUND_FIELD = 4,
    MS_ERROR_NONFINITE_INPUT = 5,
    MS_ERROR_INPUT_TOO_SHORT = 6,
    MS_ERROR_OUTPUT_TOO_SMALL = 7,
    MS_ERROR_SIZE_OVERFLOW = 8,
    MS_ERROR_ALLOCATION = 9,
    MS_ERROR_NUMERICAL = 10
};

typedef uint32_t ms_sound_field;
enum ms_sound_field_value {
    MS_SOUND_FIELD_FREE = 0,
    MS_SOUND_FIELD_DIFFUSE = 1
};

/**
 * Callers retain ownership of every caller-supplied buffer. Inputs are
 * read-only. Input and output ranges, and distinct writable output ranges,
 * must not overlap. Pointers returned by metadata accessors below are
 * immutable library-owned storage and must not be freed. Except where a
 * function explicitly documents otherwise, output values are valid only
 * when the function returns MS_OK and are unspecified after an error.
 */

/** Return the stable ABI version implemented by this library. */
MS_API uint32_t MS_CALL ms_abi_version(void);

/** Return the semantic library version as a static UTF-8 string. */
MS_API const char *MS_CALL ms_version_string(void);

/** Return a static English description of a status code. */
MS_API const char *MS_CALL ms_status_string(ms_status status);

/** Return immutable library-owned storage valid for the library lifetime. */
MS_API const double *MS_CALL ms_third_octave_centres_hz(void);

/** Return immutable library-owned storage valid for the library lifetime. */
MS_API const double *MS_CALL ms_bark_axis(void);

/**
 * Compute ISO 532-1:2017-targeted stationary loudness from 28 third-octave
 * levels (25 Hz through 12.5 kHz, dB SPL re 20 uPa).
 *
 * The first 11 values (25 Hz through 250 Hz) must not exceed 120 dB SPL;
 * higher values would extrapolate beyond the implemented low-band model data.
 *
 * specific_sone_per_bark may be NULL only when specific_capacity is zero.
 * When supplied, its capacity must be at least MS_BARK_BANDS.
 */
MS_API ms_status MS_CALL ms_loudness_from_levels(
    const double *levels_db,
    size_t level_count,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity);

/**
 * Compute stationary loudness from mono pressure samples in pascals.
 * The native signal kernel deliberately accepts exactly 48 kHz only.
 */
MS_API ms_status MS_CALL ms_loudness_stationary(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity);

/**
 * Query the number of 2 ms frames produced by time-varying loudness.
 * Unaddressable sample counts return MS_ERROR_SIZE_OVERFLOW.
 */
MS_API ms_status MS_CALL ms_loudness_time_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);

/**
 * Compute time-varying loudness. Arrays are frame-major; specific output
 * requires frame_count * MS_BARK_BANDS elements. frame_count_written is
 * populated on both success and MS_ERROR_OUTPUT_TOO_SMALL.
 */
MS_API ms_status MS_CALL ms_loudness_time(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *loudness_sone,
    size_t loudness_capacity,
    double *specific_sone_per_bark,
    size_t specific_capacity,
    size_t *frame_count_written);

/**
 * Query the 200 ms / 50% overlap frame count used for roughness.
 * Sample counts exceeding INT_MAX or addressable memory return MS_ERROR_SIZE_OVERFLOW.
 */
MS_API ms_status MS_CALL ms_roughness_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);

/**
 * Compute the experimental Daniel-Weber-targeted tabulated roughness model in
 * asper at 48 kHz. See the documented validation limits before use.
 * frame_count_written receives the required count on MS_ERROR_OUTPUT_TOO_SMALL.
 */
MS_API ms_status MS_CALL ms_roughness_dw(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    double *roughness_asper,
    size_t roughness_capacity,
    size_t *frame_count_written);

/**
 * Query the 250 ms / 50% overlap frame count used for Aures tonality.
 * Sample counts exceeding INT_MAX or addressable memory return MS_ERROR_SIZE_OVERFLOW.
 */
MS_API ms_status MS_CALL ms_tonality_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);

/**
 * Compute experimental Aures-1985-targeted tonality at 48 kHz. See the
 * documented validation limits before use.
 * frame_count_written receives the required count on MS_ERROR_OUTPUT_TOO_SMALL.
 */
MS_API ms_status MS_CALL ms_tonality_aures(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *tonality,
    size_t tonality_capacity,
    size_t *frame_count_written);

/**
 * Compute DIN 45692 sharpness from one specific-loudness pattern.
 * The input must contain exactly MS_BARK_BANDS values in sone/Bark.
 */
MS_API ms_status MS_CALL ms_sharpness_din(
    const double *specific_sone_per_bark,
    size_t specific_count,
    double *sharpness_acum);

#ifdef __cplusplus
}
#endif

#endif
