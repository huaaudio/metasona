// SPDX-License-Identifier: Apache-2.0 AND MIT
// Adapted roughness model; SQAT: copyright 1999-2023 Dik Hermes.
// Modified for the MetaSona C API; see NOTICE and THIRD_PARTY.md.

#include "ms_internal.h"
#include "ms_roughness_data.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Rewrite of the supplied Daniel--Weber reference stages. All scratch storage
 * belongs to one call. Model data provenance is recorded separately. */
typedef struct ms_roughness_scratch {
    ms_complex *spectrum;
    ms_complex *work;
    double *window;
    double *modulation;
    double *spectral_bark;
    double *spectral_audibility;
    double *upper_spreading_slope;
    double *spectral_level;
    size_t fft_size;
    double window_sum;
} ms_roughness_scratch;

static double ms_interpolate_table(const double table[][2], size_t count, double x)
{
    size_t index;
    if (x <= table[0][0]) {
        return table[0][1];
    }
    for (index = 1u; index < count; ++index) {
        if (x <= table[index][0]) {
            const double ratio = (x - table[index - 1u][0])
                / (table[index][0] - table[index - 1u][0]);
            return table[index - 1u][1]
                + ratio * (table[index][1] - table[index - 1u][1]);
        }
    }
    return table[count - 1u][1];
}

static double ms_roughness_bark(double frequency)
{
    /* Interleaved band edges and centres are already ordered by frequency. */
    double previous_frequency = 0.0;
    double previous_bark = 0.0;
    size_t point;
    for (point = 1u; point < 50u; ++point) {
        const size_t row = point / 2u;
        const int centre = (int)(point % 2u);
        const double next_frequency = MS_ROUGHNESS_BARK[row][centre ? 2 : 1];
        const double next_bark = (double)row + (centre ? 0.5 : 0.0);
        if (frequency <= next_frequency) {
            return previous_bark + (next_bark - previous_bark)
                * (frequency - previous_frequency) / (next_frequency - previous_frequency);
        }
        previous_frequency = next_frequency;
        previous_bark = next_bark;
    }
    return 24.5;
}

static double ms_modulation_weight(size_t channel, double frequency)
{
    if (frequency < 10.0) {
        return 0.0;
    }
    if (channel < 4u) {
        return frequency <= 355.0
            ? ms_interpolate_table(MS_ROUGHNESS_H2, 15u, frequency) : 0.0;
    }
    /* Preserve the supplied reference's 502 Hz cutoff (500 Hz at 5 Hz bins),
     * including for the H16/H21/H42 tables which extend to 645 Hz. */
    if (frequency > 500.0) {
        return 0.0;
    }
    if (channel < 15u) {
        return ms_interpolate_table(MS_ROUGHNESS_H5, 14u, frequency);
    }
    if (channel < 20u) {
        return ms_interpolate_table(MS_ROUGHNESS_H16, 20u, frequency);
    }
    if (channel < 41u) {
        return ms_interpolate_table(MS_ROUGHNESS_H21, 18u, frequency);
    }
    return ms_interpolate_table(MS_ROUGHNESS_H42, 19u, frequency);
}

static double ms_channel_weight(size_t channel)
{
    const double bark = 0.5 * (double)(channel + 1u);
    const size_t low = (size_t)bark;
    return MS_ROUGHNESS_CHANNEL_WEIGHT[low] + (bark - (double)low)
        * (MS_ROUGHNESS_CHANNEL_WEIGHT[low + 1u] - MS_ROUGHNESS_CHANNEL_WEIGHT[low]);
}

static double ms_roughness_correlation(const double *left, const double *right, size_t count)
{
    double cross = 0.0, left_square = 0.0, right_square = 0.0;
    double left_mean = 0.0, right_mean = 0.0;
    size_t index;
    for (index = 0u; index < count; ++index) {
        left_mean += left[index] / (double)count;
        right_mean += right[index] / (double)count;
    }
    for (index = 0u; index < count; ++index) {
        const double a = left[index] - left_mean;
        const double b = right[index] - right_mean;
        cross += a * b;
        left_square += a * a;
        right_square += b * b;
    }
    if (left_square <= 0.0 || right_square <= 0.0) {
        return 0.0;
    }
    /* Divide separately to avoid overflowing a product of squared norms. */
    return fmax(-1.0, fmin(1.0, (cross / sqrt(left_square)) / sqrt(right_square)));
}

static ms_status ms_roughness_prepare_spectrum(
    const double *frame, ms_roughness_scratch *scratch, int *audible)
{
    const size_t count = scratch->fft_size;
    const double calibration = 2.0 * pow(10.0, 91.2 / 20.0) / scratch->window_sum;
    size_t index;
    ms_status status;
    *audible = 0;
    for (index = 0u; index < count; ++index) {
        scratch->spectrum[index].re = frame[index] * scratch->window[index] * calibration;
        scratch->spectrum[index].im = 0.0;
    }
    status = ms_fft_forward(scratch->spectrum, count);
    if (status != MS_OK) {
        return status;
    }
    memset(scratch->spectral_audibility, 0, (count / 2u + 1u) * sizeof(double));
    /* True zero-based DFT bins: bin 4 is 20 Hz. The legacy port mixed
     * one-based spectrum offsets and Bark lookups. No offset is added here. */
    for (index = 4u; index <= 4000u; ++index) {
        const double frequency = 5.0 * (double)index;
        const double bark = ms_roughness_bark(frequency);
        const double ear = pow(10.0,
            ms_interpolate_table(MS_ROUGHNESS_EAR, 22u, bark) / 20.0);
        double magnitude;
        scratch->spectrum[index].re *= ear;
        scratch->spectrum[index].im *= ear;
        magnitude = hypot(scratch->spectrum[index].re, scratch->spectrum[index].im);
        scratch->spectral_bark[index] = bark;
        if (!ms_is_finite(magnitude)) {
            return MS_ERROR_NUMERICAL;
        }
        if (magnitude > 0.0) {
            const double level = 20.0 * log10(magnitude);
            scratch->spectral_level[index] = level;
            if (level > ms_interpolate_table(MS_ROUGHNESS_THRESHOLD, 27u, bark)) {
                scratch->spectral_audibility[index] = 1.0;
                scratch->upper_spreading_slope[index] = fmin(
                    -24.0 - 230.0 / frequency + 0.2 * level, 0.0);
                *audible = 1;
            }
        }
    }
    return MS_OK;
}

static ms_status ms_roughness_channel_modulation(
    size_t channel, ms_roughness_scratch *scratch, double *depth)
{
    const size_t count = scratch->fft_size;
    const double centre = 0.5 * (double)(channel + 1u);
    double mean = 0.0, square = 0.0;
    size_t index;
    ms_status status;
    memset(scratch->work, 0, count * sizeof(*scratch->work));
    for (index = 4u; index <= 4000u; ++index) {
        const double source = scratch->spectral_bark[index];
        double spreading = 0.0;
        double target = centre;
        double gain;
        if (scratch->spectral_audibility[index] == 0.0) {
            continue;
        }
        if (source > centre + 0.5) {
            target = centre + 0.5;
            spreading = -27.0 * (source - target);
        } else if (source < centre - 0.5) {
            target = centre - 0.5;
            spreading = scratch->upper_spreading_slope[index] * (target - source);
        }
        if (spreading < 0.0 && scratch->spectral_level[index] + spreading
            <= ms_interpolate_table(MS_ROUGHNESS_THRESHOLD, 27u, target)) {
            continue;
        }
        gain = pow(10.0, spreading / 20.0);
        scratch->work[index].re = gain * scratch->spectrum[index].re;
        scratch->work[index].im = gain * scratch->spectrum[index].im;
        scratch->work[count - index].re = scratch->work[index].re;
        scratch->work[count - index].im = -scratch->work[index].im;
    }
    status = ms_fft_inverse(scratch->work, count);
    if (status != MS_OK) {
        return status;
    }
    for (index = 0u; index < count; ++index) {
        mean += fabs(scratch->work[index].re) / (double)count;
    }
    if (!ms_is_finite(mean)) {
        return MS_ERROR_NUMERICAL;
    }
    *depth = 0.0;
    if (mean <= 0.0) {
        memset(scratch->modulation + channel * count, 0, count * sizeof(double));
        return MS_OK;
    }
    /* Normalizing before the envelope transform makes depth and correlation
     * scale invariant and avoids large intermediate envelope energies. */
    for (index = 0u; index < count; ++index) {
        scratch->work[index].re = fabs(scratch->work[index].re) / mean - 1.0;
        scratch->work[index].im = 0.0;
    }
    status = ms_fft_forward(scratch->work, count);
    if (status != MS_OK) {
        return status;
    }
    for (index = 0u; index <= count / 2u; ++index) {
        const double weight = ms_modulation_weight(channel, (double)index * 5.0);
        scratch->work[index].re *= weight;
        scratch->work[index].im *= weight;
        if (index > 0u && index < count / 2u) {
            scratch->work[count - index].re *= weight;
            scratch->work[count - index].im *= weight;
        }
    }
    status = ms_fft_inverse(scratch->work, count);
    if (status != MS_OK) {
        return status;
    }
    for (index = 0u; index < count; ++index) {
        const double value = scratch->work[index].re;
        scratch->modulation[channel * count + index] = value;
        square += value * value / (double)count;
    }
    if (!ms_is_finite(square)) {
        return MS_ERROR_NUMERICAL;
    }
    *depth = fmin(1.0, sqrt(square));
    return MS_OK;
}

static ms_status ms_roughness_one_frame(
    const double *frame, ms_roughness_scratch *scratch, double *roughness)
{
    double depth[47], correlation[45];
    double sum = 0.0;
    size_t channel;
    int audible;
    ms_status status = ms_roughness_prepare_spectrum(frame, scratch, &audible);
    if (status != MS_OK) {
        return status;
    }
    if (!audible) {
        *roughness = 0.0;
        return MS_OK;
    }
    for (channel = 0u; channel < 47u; ++channel) {
        status = ms_roughness_channel_modulation(channel, scratch, &depth[channel]);
        if (status != MS_OK) {
            return status;
        }
    }
    for (channel = 0u; channel < 45u; ++channel) {
        correlation[channel] = ms_roughness_correlation(
            scratch->modulation + channel * scratch->fft_size,
            scratch->modulation + (channel + 2u) * scratch->fft_size, scratch->fft_size);
    }
    for (channel = 0u; channel < 47u; ++channel) {
        double term = depth[channel];
        if (channel < 45u) {
            term *= correlation[channel];
        }
        if (channel >= 2u) {
            term *= correlation[channel - 2u];
        }
        /* The empirical channel weight is outside the square. */
        sum += ms_channel_weight(channel) * term * term;
    }
    *roughness = 0.25 * sum; /* CAL=0.5, Bark spacing=0.5; no fitted factor. */
    return ms_is_finite(*roughness) ? MS_OK : MS_ERROR_NUMERICAL;
}

ms_status ms_roughness_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count)
{
    const size_t frame_size = 9600u;
    const size_t hop_size = 4800u;
    size_t result;
    if (frame_count == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (sample_rate_hz != MS_SAMPLE_RATE_HZ) {
        return MS_ERROR_INVALID_SAMPLE_RATE;
    }
    if (sample_count > (size_t)PTRDIFF_MAX / sizeof(double)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    if (sample_count < frame_size) {
        return MS_ERROR_INPUT_TOO_SHORT;
    }
    result = 1u + (sample_count - frame_size) / hop_size;
    *frame_count = result;
    return MS_OK;
}

ms_status ms_roughness_dw(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    double *roughness_asper,
    size_t roughness_capacity,
    size_t *frame_count_written)
{
    const size_t frame_size = 9600u;
    const size_t hop_size = 4800u;
    const size_t channel_count = 47u;
    ms_roughness_scratch scratch;
    double *temporary_results = NULL;
    double *spectral_storage = NULL;
    size_t frame_count;
    const size_t fft_size = frame_size;
    size_t complex_bytes;
    size_t result_bytes;
    size_t modulation_count;
    size_t modulation_bytes;
    size_t spectral_count;
    size_t spectral_bytes;
    size_t frame;
    size_t index;
    ms_status status;

    if (pressure_pa == NULL || roughness_asper == NULL
        || frame_count_written == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    status = ms_roughness_frame_count(
        sample_count, sample_rate_hz, &frame_count);
    if (status != MS_OK) {
        return status;
    }
    if (roughness_capacity < frame_count) {
        *frame_count_written = frame_count;
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    if (!ms_checked_multiply(fft_size, sizeof(ms_complex), &complex_bytes)
        || !ms_checked_multiply(frame_count, sizeof(double), &result_bytes)
        || !ms_checked_multiply(channel_count, frame_size, &modulation_count)
        || !ms_checked_multiply(modulation_count, sizeof(double), &modulation_bytes)
        || !ms_checked_multiply(fft_size / 2u + 1u, 4u, &spectral_count)
        || !ms_checked_multiply(spectral_count, sizeof(double), &spectral_bytes)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    status = ms_validate_signal(pressure_pa, sample_count, sample_rate_hz);
    if (status != MS_OK) {
        return status;
    }

    memset(&scratch, 0, sizeof(scratch));
    scratch.fft_size = fft_size;
    scratch.spectrum = (ms_complex *)calloc(fft_size, sizeof(*scratch.spectrum));
    scratch.work = (ms_complex *)calloc(fft_size, sizeof(*scratch.work));
    scratch.window = (double *)calloc(frame_size, sizeof(*scratch.window));
    scratch.modulation = (double *)calloc(
        modulation_count, sizeof(*scratch.modulation));
    spectral_storage = (double *)calloc(
        spectral_count, sizeof(*spectral_storage));
    temporary_results = (double *)calloc(
        frame_count, sizeof(*temporary_results));
    if (scratch.spectrum == NULL || scratch.work == NULL
        || scratch.window == NULL || scratch.modulation == NULL
        || spectral_storage == NULL || temporary_results == NULL) {
        free(scratch.spectrum);
        free(scratch.work);
        free(scratch.window);
        free(scratch.modulation);
        free(spectral_storage);
        free(temporary_results);
        return MS_ERROR_ALLOCATION;
    }
    scratch.spectral_bark = spectral_storage;
    scratch.spectral_audibility = spectral_storage + fft_size / 2u + 1u;
    scratch.upper_spreading_slope = scratch.spectral_audibility
                                   + fft_size / 2u + 1u;
    scratch.spectral_level = scratch.upper_spreading_slope + fft_size / 2u + 1u;
    for (index = 0u; index < frame_size; ++index) {
        const double phase = 2.0 * MS_PI * (double)index / (double)(frame_size - 1u);
        const double window = 0.42 - 0.5 * cos(phase) + 0.08 * cos(2.0 * phase);
        scratch.window[index] = window;
        scratch.window_sum += window;
    }

    status = MS_OK;
    for (frame = 0u; frame < frame_count; ++frame) {
        status = ms_roughness_one_frame(
            pressure_pa + frame * hop_size,
            &scratch,
            &temporary_results[frame]);
        if (status != MS_OK) {
            break;
        }
    }
    if (status == MS_OK) {
        memcpy(roughness_asper, temporary_results, result_bytes);
        *frame_count_written = frame_count;
    }
    free(scratch.spectrum);
    free(scratch.work);
    free(scratch.window);
    free(scratch.modulation);
    free(spectral_storage);
    free(temporary_results);
    return status;
}
