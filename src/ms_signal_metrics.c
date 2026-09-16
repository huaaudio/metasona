// SPDX-License-Identifier: Apache-2.0
// MetaSona authors: Jiahua Zhang and Codex, September 2026.
// MetaSona API boundary for the reused psychohelperc kernels.
#include "ms_internal.h"
#include "roughness_dw.h"
#include "tonality_aures1985.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

static ms_status frame_count_for(size_t samples, uint32_t rate,
    size_t window, size_t hop, size_t *count)
{
    if (count == NULL) return MS_ERROR_NULL_POINTER;
    if (rate != MS_SAMPLE_RATE_HZ) return MS_ERROR_INVALID_SAMPLE_RATE;
    if (samples > (size_t)INT_MAX || samples > (size_t)PTRDIFF_MAX / sizeof(double))
        return MS_ERROR_SIZE_OVERFLOW;
    if (samples < window) return MS_ERROR_INPUT_TOO_SHORT;
    *count = 1u + (samples - window) / hop;
    return MS_OK;
}

ms_status ms_roughness_frame_count(size_t samples, uint32_t rate, size_t *count)
{
    return frame_count_for(samples, rate, 9600u, 4800u, count);
}

ms_status ms_tonality_frame_count(size_t samples, uint32_t rate, size_t *count)
{
    return frame_count_for(samples, rate, 12000u, 6000u, count);
}

static ms_status calculate(const double *pressure, size_t samples, uint32_t rate,
    ms_sound_field field, int tonal, double *output, size_t capacity, size_t *written)
{
    size_t count, bytes, i;
    int kernel_count, result;
    double *temporary;
    ms_status status;
    if (pressure == NULL || output == NULL || written == NULL)
        return MS_ERROR_NULL_POINTER;
    status = tonal ? ms_tonality_frame_count(samples, rate, &count)
                   : ms_roughness_frame_count(samples, rate, &count);
    if (status != MS_OK) return status;
    if (tonal && !ms_is_valid_field(field)) return MS_ERROR_INVALID_SOUND_FIELD;
    if (capacity < count) {
        *written = count;
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    status = ms_validate_signal(pressure, samples, rate);
    if (status != MS_OK) return status;
    if (!ms_checked_multiply(count, sizeof(double), &bytes)) return MS_ERROR_SIZE_OVERFLOW;
    temporary = (double *)malloc(bytes);
    if (temporary == NULL) return MS_ERROR_ALLOCATION;
    kernel_count = (int)count;
    result = tonal
        ? tonality_aures1985(pressure, (int)samples, (double)rate, (int)field, 0.0,
            temporary, NULL, NULL, NULL, &kernel_count, NULL)
        : roughness_dw(pressure, (int)samples, (double)rate,
            temporary, NULL, &kernel_count);
    if (result == RoughnessErrorMemoryAlloc
        || (tonal && result == TonalityAures1985ErrorFFTPlan))
        status = MS_ERROR_ALLOCATION;
    else if (result < 0 || (size_t)result != count || (size_t)kernel_count != count)
        status = MS_ERROR_NUMERICAL;
    else {
        for (i = 0u; i < count; ++i) {
            if (!ms_is_finite(temporary[i]) || temporary[i] < 0.0) {
                status = MS_ERROR_NUMERICAL;
                break;
            }
        }
    }
    if (status == MS_OK) {
        memcpy(output, temporary, bytes);
        *written = count;
    }
    free(temporary);
    return status;
}

ms_status ms_roughness_dw(const double *pressure, size_t samples, uint32_t rate,
    double *output, size_t capacity, size_t *written)
{
    return calculate(pressure, samples, rate, MS_SOUND_FIELD_FREE, 0, output, capacity, written);
}

ms_status ms_tonality_aures(const double *pressure, size_t samples, uint32_t rate,
    ms_sound_field field, double *output, size_t capacity, size_t *written)
{
    return calculate(pressure, samples, rate, field, 1, output, capacity, written);
}
