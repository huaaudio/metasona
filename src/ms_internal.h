// SPDX-License-Identifier: Apache-2.0

#ifndef METASONA_MS_INTERNAL_H
#define METASONA_MS_INTERNAL_H

#include "metasona/metasona.h"

#include <stddef.h>
#include <stdint.h>

#define MS_REFERENCE_PRESSURE_PA 0.00002
#define MS_PI 3.14159265358979323846264338327950288

typedef struct ms_complex {
    double re;
    double im;
} ms_complex;

int ms_is_valid_field(ms_sound_field field);
int ms_is_finite(double value);
ms_status ms_validate_signal(const double *samples, size_t count, uint32_t rate);
int ms_checked_multiply(size_t left, size_t right, size_t *result);
int ms_next_power_of_two(size_t value, size_t *result);
ms_status ms_fft_forward(ms_complex *values, size_t count);
ms_status ms_fft_inverse(ms_complex *values, size_t count);
double ms_hz_to_bark(double frequency_hz);
double ms_bark_to_hz(double bark);

ms_status ms_levels_from_signal_48k(
    const double *samples,
    size_t count,
    double levels_db[MS_THIRD_OCTAVE_BANDS]);

#endif
