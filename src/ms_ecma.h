// SPDX-License-Identifier: GPL-3.0-only
// MetaSona author: Jiahua Zhang, 2026.
#ifndef METASONA_MS_ECMA_H
#define METASONA_MS_ECMA_H

#include "ms_internal.h"

#define MS_ECMA_BANDS 53u

/* Shared Section 5 primitives. The edition selects 2022 rounded coefficients
 * (MoSQITo comparison) or 2025 coefficients/calibration (SQAT/ECMA target).
 * These are internal building blocks, not the Section 8 loudness metric. */
double ms_ecma_centre(size_t band);
double ms_ecma_basis(double rms_pa, size_t band, unsigned edition);
ms_status ms_ecma_prepare(const double *input, size_t count, size_t block,
    size_t hop, int pad_end, unsigned edition, ms_sound_field field,
    double **output, size_t *output_count);
ms_status ms_ecma_band(const double *input, size_t count, size_t band,
    double *output);
ms_status ms_ecma_basis_frames(const double *bandpass, size_t count,
    size_t start, size_t block, size_t hop, size_t frames, size_t band,
    unsigned edition, int mosqito_indexing, double *output);

ms_status ms_ecma_tonal_count(size_t count, size_t *frames);
/* Five frame-major planes: tonal loudness, noise loudness, tonal frequency,
 * specific tonality, combined specific loudness. Each is frames * 53. */
ms_status ms_ecma_tonal(const double *input, size_t count, ms_sound_field field,
    double *out, size_t capacity, size_t *written);

ms_status ms_ecma_rough_count(size_t count, size_t *frames);
/* Frame-major specific roughness: frames * 53, 50 Hz output grid. */
ms_status ms_ecma_rough(const double *input, size_t count, ms_sound_field field,
    double *out, size_t capacity, size_t *written);

#endif
