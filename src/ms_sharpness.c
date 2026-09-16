// SPDX-License-Identifier: Apache-2.0
// MetaSona author: Jiahua Zhang, September 2026.
#include "ms_internal.h"
#include <math.h>

ms_status ms_sharpness_din(
    const double *specific_sone_per_bark,
    size_t specific_count,
    double *sharpness_acum)
{
    double maximum = 0.0;
    double denominator = 0.0;
    double numerator = 0.0;
    double result;
    size_t index;
    if (specific_sone_per_bark == NULL || sharpness_acum == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (specific_count != MS_BARK_BANDS) {
        return MS_ERROR_INVALID_ARGUMENT;
    }
    for (index = 0u; index < MS_BARK_BANDS; ++index) {
        const double value = specific_sone_per_bark[index];
        if (!ms_is_finite(value)) {
            return MS_ERROR_NONFINITE_INPUT;
        }
        if (value < 0.0) {
            return MS_ERROR_INVALID_ARGUMENT;
        }
        if (value > maximum) {
            maximum = value;
        }
    }
    if (maximum == 0.0) {
        *sharpness_acum = 0.0;
        return MS_OK;
    }
    for (index = 0u; index < MS_BARK_BANDS; ++index) {
        const double z = 0.1 * (double)(index + 1u);
        /* Common scale cancels in the ratio; prevents overflow without
         * changing the loudness pattern or calibrating the result. */
        const double scaled_loudness = specific_sone_per_bark[index] / maximum;
        const double weighting = z <= 15.8
            ? 1.0 : 0.15 * exp(0.42 * (z - 15.8)) + 0.85;
        denominator += scaled_loudness;
        numerator += scaled_loudness * weighting * z;
    }
    if (denominator <= 0.0) {
        return MS_ERROR_NUMERICAL;
    }
    result = 0.11 * numerator / denominator;
    if (!ms_is_finite(result) || result < 0.0) {
        return MS_ERROR_NUMERICAL;
    }
    *sharpness_acum = result;
    return MS_OK;
}
