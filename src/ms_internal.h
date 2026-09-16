// SPDX-License-Identifier: Apache-2.0
// MetaSona author: Jiahua Zhang, September 2026.

#ifndef METASONA_MS_INTERNAL_H
#define METASONA_MS_INTERNAL_H

#include "metasona/metasona.h"

#include <stddef.h>
#include <stdint.h>

#define MS_REFERENCE_PRESSURE_PA 0.00002
#define MS_PI 3.14159265358979323846264338327950288

int ms_is_valid_field(ms_sound_field field);
int ms_is_finite(double value);
ms_status ms_validate_signal(const double *samples, size_t count, uint32_t rate);
int ms_checked_multiply(size_t left, size_t right, size_t *result);
ms_status ms_levels_from_signal_48k(
    const double *samples,
    size_t count,
    double levels_db[MS_THIRD_OCTAVE_BANDS]);

ms_status ms_levels_from_signal_skip_48k(
    const double *samples, size_t count, size_t skip,
    double levels_db[MS_THIRD_OCTAVE_BANDS]);

#endif
