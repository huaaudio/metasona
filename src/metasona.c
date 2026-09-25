// SPDX-License-Identifier: Apache-2.0
// MetaSona author: Jiahua Zhang, September 2026.

#include "ms_internal.h"

#include <math.h>
#include <float.h>
#include <stdint.h>

static const double MS_CENTRES_HZ[MS_THIRD_OCTAVE_BANDS] = {
    25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0,
    125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0,
    630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0,
    3150.0, 4000.0, 5000.0, 6300.0, 8000.0, 10000.0, 12500.0
};

static const double MS_BARK_VALUES[MS_BARK_BANDS] = {
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
    1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
    2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0,
    3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0,
    4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0,
    5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0,
    6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0,
    7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0,
    8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0,
    9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0,
    10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0,
    11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0,
    12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0,
    13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0,
    14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0,
    15.1, 15.2, 15.3, 15.4, 15.5, 15.6, 15.7, 15.8, 15.9, 16.0,
    16.1, 16.2, 16.3, 16.4, 16.5, 16.6, 16.7, 16.8, 16.9, 17.0,
    17.1, 17.2, 17.3, 17.4, 17.5, 17.6, 17.7, 17.8, 17.9, 18.0,
    18.1, 18.2, 18.3, 18.4, 18.5, 18.6, 18.7, 18.8, 18.9, 19.0,
    19.1, 19.2, 19.3, 19.4, 19.5, 19.6, 19.7, 19.8, 19.9, 20.0,
    20.1, 20.2, 20.3, 20.4, 20.5, 20.6, 20.7, 20.8, 20.9, 21.0,
    21.1, 21.2, 21.3, 21.4, 21.5, 21.6, 21.7, 21.8, 21.9, 22.0,
    22.1, 22.2, 22.3, 22.4, 22.5, 22.6, 22.7, 22.8, 22.9, 23.0,
    23.1, 23.2, 23.3, 23.4, 23.5, 23.6, 23.7, 23.8, 23.9, 24.0
};

uint32_t ms_abi_version(void)
{
    return MS_ABI_VERSION;
}

const char *ms_version_string(void)
{
    return "0.2.2";
}

const char *ms_status_string(ms_status status)
{
    switch (status) {
    case MS_OK: return "success";
    case MS_ERROR_NULL_POINTER: return "null pointer";
    case MS_ERROR_INVALID_ARGUMENT: return "invalid argument";
    case MS_ERROR_INVALID_SAMPLE_RATE: return "invalid sample rate";
    case MS_ERROR_INVALID_SOUND_FIELD: return "invalid sound field";
    case MS_ERROR_NONFINITE_INPUT: return "input contains NaN or infinity";
    case MS_ERROR_INPUT_TOO_SHORT: return "input is too short";
    case MS_ERROR_OUTPUT_TOO_SMALL: return "output buffer is too small";
    case MS_ERROR_SIZE_OVERFLOW: return "size calculation overflow";
    case MS_ERROR_ALLOCATION: return "memory allocation failed";
    case MS_ERROR_NUMERICAL: return "numerical calculation failed";
    default: return "unknown status";
    }
}

const double *ms_third_octave_centres_hz(void)
{
    return MS_CENTRES_HZ;
}

const double *ms_bark_axis(void)
{
    return MS_BARK_VALUES;
}

int ms_is_valid_field(ms_sound_field field)
{
    return field == MS_SOUND_FIELD_FREE || field == MS_SOUND_FIELD_DIFFUSE;
}

int ms_is_finite(double value)
{
    return value <= DBL_MAX && value >= -DBL_MAX;
}

ms_status ms_validate_signal(const double *samples, size_t count, uint32_t rate)
{
    size_t index;
    if (samples == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (count == 0u) {
        return MS_ERROR_INPUT_TOO_SHORT;
    }
    if (rate != MS_SAMPLE_RATE_HZ) {
        return MS_ERROR_INVALID_SAMPLE_RATE;
    }
    if (count > (size_t)PTRDIFF_MAX / sizeof(*samples)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    for (index = 0; index < count; ++index) {
        if (!ms_is_finite(samples[index])) {
            return MS_ERROR_NONFINITE_INPUT;
        }
    }
    return MS_OK;
}

int ms_checked_multiply(size_t left, size_t right, size_t *result)
{
    if (result == NULL) {
        return 0;
    }
    if (right != 0u && left > SIZE_MAX / right) {
        return 0;
    }
    *result = left * right;
    return 1;
}
