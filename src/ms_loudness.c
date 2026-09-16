// SPDX-License-Identifier: Apache-2.0 AND BSD-3-Clause
// MetaSona adaptation: Jiahua Zhang and Codex, September 2026.
// Adapted loudness processing; SQAT: Copyright (c) <2015>, <Ella Manor>.
// Modified for the MetaSona C API; see NOTICE and THIRD_PARTY.md.

#include "ms_internal.h"
#include "ms_loudness_data.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Zwicker model tables; upstream sources and modifications: THIRD_PARTY.md. */
static const double MS_RANGE_LIMITS[8] = {
    45.0, 55.0, 65.0, 71.0, 80.0, 90.0, 100.0, 120.0
};

static const double MS_LOW_FREQUENCY_CORRECTION[8][11] = {
    {-32.0, -24.0, -16.0, -10.0, -5.0, 0.0, -7.0, -3.0, 0.0, -2.0, 0.0},
    {-29.0, -22.0, -15.0, -10.0, -4.0, 0.0, -7.0, -2.0, 0.0, -2.0, 0.0},
    {-27.0, -19.0, -14.0,  -9.0, -4.0, 0.0, -6.0, -2.0, 0.0, -2.0, 0.0},
    {-25.0, -17.0, -12.0,  -9.0, -3.0, 0.0, -5.0, -2.0, 0.0, -2.0, 0.0},
    {-23.0, -16.0, -11.0,  -7.0, -3.0, 0.0, -4.0, -1.0, 0.0, -1.0, 0.0},
    {-20.0, -14.0, -10.0,  -6.0, -3.0, 0.0, -4.0, -1.0, 0.0, -1.0, 0.0},
    {-18.0, -12.0,  -9.0,  -6.0, -2.0, 0.0, -3.0, -1.0, 0.0, -1.0, 0.0},
    {-15.0, -10.0,  -8.0,  -4.0, -2.0, 0.0, -3.0, -1.0, 0.0, -1.0, 0.0}
};

static const double MS_THRESHOLD[20] = {
    30.0, 18.0, 12.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 3.0,
    3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0
};

static const double MS_EAR_TRANSMISSION[20] = {
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    -0.5, -1.6, -3.2, -5.4, -5.6, -4.0, -1.5, 2.0, 5.0, 12.0
};

static const double MS_DIFFUSE_CORRECTION[20] = {
    0.0, 0.0, 0.5, 0.9, 1.2, 1.6, 2.3, 2.8, 3.0, 2.0,
    0.0, -1.4, -2.0, -1.9, -1.0, 0.5, 3.0, 4.0, 4.3, 4.0
};

static const double MS_CRITICAL_BAND_CORRECTION[20] = {
    -0.25, -0.6, -0.8, -0.8, -0.5, 0.0, 0.5, 1.1, 1.5, 1.7,
    1.8, 1.8, 1.7, 1.6, 1.4, 1.2, 0.8, 0.5, 0.0, -0.5
};

static const double MS_BAND_UPPER_BARK[21] = {
    0.9, 1.8, 2.8, 3.5, 4.4, 5.4, 6.6, 7.9, 9.2, 10.6, 12.3,
    13.8, 15.2, 16.7, 18.1, 19.3, 20.6, 21.8, 22.7, 23.6, 24.0
};

static const double MS_SLOPE_RANGES[18] = {
    21.5, 18.0, 15.1, 11.5, 9.0, 6.1, 4.4, 3.1, 2.13,
    1.36, 0.82, 0.42, 0.30, 0.22, 0.15, 0.10, 0.035, 0.0
};

/* Row 12, column 4 is 0.22 in supplied ISO_532-1.c (not MoSQITo's 0.24).
 * Retained from that source, not adjusted to make a fixture pass. */
static const double MS_UPPER_SLOPES[18][8] = {
    {13.0, 8.2, 6.3, 5.5, 5.5, 5.5, 5.5, 5.5},
    {9.0, 7.5, 6.0, 5.1, 4.5, 4.5, 4.5, 4.5},
    {7.8, 6.7, 5.6, 4.9, 4.4, 3.9, 3.9, 3.9},
    {6.2, 5.4, 4.6, 4.0, 3.5, 3.2, 3.2, 3.2},
    {4.5, 3.8, 3.6, 3.2, 2.9, 2.7, 2.7, 2.7},
    {3.7, 3.0, 2.8, 2.35, 2.2, 2.2, 2.2, 2.2},
    {2.9, 2.3, 2.1, 1.9, 1.8, 1.7, 1.7, 1.7},
    {2.4, 1.7, 1.5, 1.35, 1.3, 1.3, 1.3, 1.3},
    {1.95, 1.45, 1.3, 1.15, 1.1, 1.1, 1.1, 1.1},
    {1.5, 1.2, 0.94, 0.86, 0.82, 0.82, 0.82, 0.82},
    {0.72, 0.67, 0.64, 0.63, 0.62, 0.62, 0.62, 0.62},
    {0.59, 0.53, 0.51, 0.50, 0.42, 0.42, 0.42, 0.42},
    {0.40, 0.33, 0.26, 0.24, 0.22, 0.22, 0.22, 0.22},
    {0.27, 0.21, 0.20, 0.18, 0.17, 0.17, 0.17, 0.17},
    {0.16, 0.15, 0.14, 0.12, 0.11, 0.11, 0.11, 0.11},
    {0.12, 0.11, 0.10, 0.08, 0.08, 0.08, 0.08, 0.08},
    {0.09, 0.08, 0.07, 0.06, 0.06, 0.06, 0.06, 0.05},
    {0.06, 0.05, 0.03, 0.02, 0.02, 0.02, 0.02, 0.02}
};

static size_t ms_slope_index(double specific)
{
    size_t index;
    for (index = 0u; index + 1u < 18u; ++index) {
        if (specific > MS_SLOPE_RANGES[index] + 1.0e-12) {
            return index;
        }
    }
    return 17u;
}

static void ms_core_loudness(
    const double levels_db[MS_THIRD_OCTAVE_BANDS],
    ms_sound_field field,
    double core[21])
{
    double low_intensity[11];
    double critical_level[3];
    size_t band;
    size_t group;

    for (band = 0u; band < 11u; ++band) {
        size_t range = 0u;
        while (range + 1u < 8u
               && levels_db[band]
                    > MS_RANGE_LIMITS[range]
                    - MS_LOW_FREQUENCY_CORRECTION[range][band]) {
            ++range;
        }
        low_intensity[band] = pow(
            10.0,
            0.1 * (levels_db[band] + MS_LOW_FREQUENCY_CORRECTION[range][band]));
    }

    for (group = 0u; group < 3u; ++group) {
        static const size_t starts[3] = {0u, 6u, 9u};
        static const size_t ends[3] = {6u, 9u, 11u};
        double intensity = 0.0;
        for (band = starts[group]; band < ends[group]; ++band) {
            intensity += low_intensity[band];
        }
        critical_level[group] = intensity > 0.0 ? 10.0 * log10(intensity) : -120.0;
    }

    for (band = 0u; band < 20u; ++band) {
        double effective_level = levels_db[band + 8u];
        if (band < 3u) {
            effective_level = critical_level[band];
        }
        effective_level -= MS_EAR_TRANSMISSION[band];
        if (field == MS_SOUND_FIELD_DIFFUSE) {
            effective_level += MS_DIFFUSE_CORRECTION[band];
        }
        core[band] = 0.0;
        if (effective_level > MS_THRESHOLD[band]) {
            const double corrected_level = effective_level - MS_CRITICAL_BAND_CORRECTION[band];
            const double threshold_factor = 0.0635 * pow(10.0, 0.025 * MS_THRESHOLD[band]);
            const double excitation = pow(
                0.75 + 0.25 * pow(10.0, 0.1 * (corrected_level - MS_THRESHOLD[band])),
                0.25) - 1.0;
            core[band] = threshold_factor * excitation;
            if (core[band] < 0.0) {
                core[band] = 0.0;
            }
        }
    }
    core[20] = 0.0;
    {
        const double correction = 0.4 + 0.32 * pow(core[0], 0.2);
        if (correction <= 1.0) {
            core[0] *= correction;
        }
    }
}

static double ms_specific_loudness(
    const double core[21],
    double specific[MS_BARK_BANDS])
{
    double z1 = 0.0;
    double n1 = 0.0;
    double total = 0.0;
    size_t output_index = 0u;
    size_t band;
    for (band = 0u; band < 21u; ++band) {
        const double limit = MS_BAND_UPPER_BARK[band];
        const double target = core[band];
        size_t slope_column = band == 0u ? 0u : band - 1u;
        if (slope_column > 7u) {
            slope_column = 7u;
        }
        while (z1 < limit - 1.0e-12) {
            double z2;
            double n2;
            if (n1 <= target) {
                n1 = target;
                z2 = limit;
                n2 = target;
            } else {
                const size_t slope_row = ms_slope_index(n1);
                const double slope = MS_UPPER_SLOPES[slope_row][slope_column];
                n2 = MS_SLOPE_RANGES[slope_row];
                if (n2 < target) {
                    n2 = target;
                }
                z2 = z1 + (n1 - n2) / slope;
                if (z2 > limit) {
                    z2 = limit;
                    n2 = n1 - (z2 - z1) * slope;
                }
            }
            total += 0.5 * (n1 + n2) * (z2 - z1);
            while (output_index < MS_BARK_BANDS) {
                const double output_z = 0.1 * (double)(output_index + 1u);
                if (output_z > z2 + 1.0e-10) {
                    break;
                }
                if (z2 > z1) {
                    const double fraction = (output_z - z1) / (z2 - z1);
                    specific[output_index] = n1 + fraction * (n2 - n1);
                } else {
                    specific[output_index] = n2;
                }
                if (specific[output_index] < 0.0) {
                    specific[output_index] = 0.0;
                }
                ++output_index;
            }
            z1 = z2;
            n1 = n2;
        }
    }
    while (output_index < MS_BARK_BANDS) {
        specific[output_index++] = 0.0;
    }
    return total;
}

ms_status ms_loudness_from_levels(
    const double *levels_db,
    size_t level_count,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity)
{
    double core[21];
    double temporary[MS_BARK_BANDS];
    double result;
    size_t index;
    if (levels_db == NULL || loudness_sone == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (level_count != MS_THIRD_OCTAVE_BANDS) {
        return MS_ERROR_INVALID_ARGUMENT;
    }
    if (!ms_is_valid_field(field)) {
        return MS_ERROR_INVALID_SOUND_FIELD;
    }
    if ((specific_sone_per_bark == NULL && specific_capacity != 0u)
        || (specific_sone_per_bark != NULL && specific_capacity < MS_BARK_BANDS)) {
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    for (index = 0u; index < MS_THIRD_OCTAVE_BANDS; ++index) {
        if (!ms_is_finite(levels_db[index])) {
            return MS_ERROR_NONFINITE_INPUT;
        }
        if (index < 11u && levels_db[index] > 120.0) {
            return MS_ERROR_INVALID_ARGUMENT;
        }
    }
    ms_core_loudness(levels_db, field, core);
    for (index = 0u; index < 21u; ++index) {
        if (!ms_is_finite(core[index])) return MS_ERROR_NUMERICAL;
    }
    result = ms_specific_loudness(core, temporary);
    if (!ms_is_finite(result)) {
        return MS_ERROR_NUMERICAL;
    }
    *loudness_sone = result;
    if (specific_sone_per_bark != NULL) {
        memcpy(specific_sone_per_bark, temporary, sizeof(temporary));
    }
    return MS_OK;
}


/* Per-call direct-form-II state. The coefficient data and processing order
 * follow the referenced model; ownership and error handling are per call. */
typedef struct ms_filter_state {
    double delay[3][2];
} ms_filter_state;

static double ms_filter_sample(ms_filter_state *state, size_t band, double input)
{
    size_t stage;
    for (stage = 0u; stage < 3u; ++stage) {
        const double *reference = MS_FILTER_REFERENCE[stage];
        const double *delta = MS_FILTER_DELTA[band][stage];
        const double first = state->delay[stage][0];
        const double second = state->delay[stage][1];
        const double next = input * MS_FILTER_GAIN[band][stage]
            - (reference[4] - delta[4]) * first
            - (reference[5] - delta[5]) * second;
        input = (reference[0] - delta[0]) * next
            + (reference[1] - delta[1]) * first
            + (reference[2] - delta[2]) * second;
        state->delay[stage][1] = first;
        state->delay[stage][0] = next;
    }
    return input;
}

static double ms_energy_level(double energy)
{
    /* Preserve the reference's absolute 1e-12 Pa^2 numerical floor. */
    return 10.0 * (log10(energy + 1.0e-12) - log10(4.0e-10));
}

ms_status ms_levels_from_signal_48k(
    const double *samples,
    size_t count,
    double levels_db[MS_THIRD_OCTAVE_BANDS])
{
    return ms_levels_from_signal_skip_48k(samples, count, 0u, levels_db);
}

ms_status ms_levels_from_signal_skip_48k(
    const double *samples, size_t count, size_t skip,
    double levels_db[MS_THIRD_OCTAVE_BANDS])
{
    double result[MS_THIRD_OCTAVE_BANDS];
    size_t band;
    if (samples == NULL || levels_db == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (count < 256u) {
        return MS_ERROR_INPUT_TOO_SHORT;
    }
    if (skip >= count) {
        return MS_ERROR_INVALID_ARGUMENT;
    }
    for (band = 0u; band < MS_THIRD_OCTAVE_BANDS; ++band) {
        ms_filter_state state = {{{0.0}}};
        double energy = 0.0;
        size_t sample;
        for (sample = 0u; sample < count; ++sample) {
            const double filtered = ms_filter_sample(&state, band, samples[sample]);
            if (sample >= skip) {
                energy += filtered * filtered;
            }
            if (!ms_is_finite(energy)) {
                return MS_ERROR_NUMERICAL;
            }
        }
        result[band] = ms_energy_level(energy / (double)(count - skip));
    }
    memcpy(levels_db, result, sizeof(result));
    return MS_OK;
}

ms_status ms_loudness_stationary(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity)
{
    double levels[MS_THIRD_OCTAVE_BANDS];
    ms_status status;
    if (loudness_sone == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (!ms_is_valid_field(field)) {
        return MS_ERROR_INVALID_SOUND_FIELD;
    }
    if ((specific_sone_per_bark == NULL && specific_capacity != 0u)
        || (specific_sone_per_bark != NULL && specific_capacity < MS_BARK_BANDS)) {
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    status = ms_validate_signal(pressure_pa, sample_count, sample_rate_hz);
    if (status == MS_OK) {
        status = ms_levels_from_signal_48k(pressure_pa, sample_count, levels);
    }
    if (status != MS_OK) {
        return status;
    }
    return ms_loudness_from_levels(levels, MS_THIRD_OCTAVE_BANDS, field,
        loudness_sone, specific_sone_per_bark, specific_capacity);
}

/* Two coupled decay states, integrated at 48 kHz between 2 kHz core values.
 * A rising core is transmitted immediately; the second state stores the
 * slower history needed to distinguish short and sustained sounds. */
typedef struct ms_decay_state {
    double output;
    double memory;
} ms_decay_state;

static void ms_decay_coefficients(double b[6])
{
    const double short_s = 0.005;
    const double long_s = 0.015;
    const double variable_s = 0.075;
    const double dt = 1.0 / 48000.0;
    const double p = (variable_s + long_s) / (variable_s * short_s);
    const double q = 1.0 / (short_s * variable_s);
    const double root = sqrt(p * p / 4.0 - q);
    const double first = -p / 2.0 + root;
    const double second = -p / 2.0 - root;
    const double denominator = variable_s * (first - second);
    const double e1 = exp(first * dt);
    const double e2 = exp(second * dt);
    b[0] = (e1 - e2) / denominator;
    b[1] = ((variable_s * second + 1.0) * e1
            - (variable_s * first + 1.0) * e2) / denominator;
    b[2] = ((variable_s * first + 1.0) * e1
            - (variable_s * second + 1.0) * e2) / denominator;
    b[3] = (variable_s * first + 1.0) * (variable_s * second + 1.0)
            * (e1 - e2) / denominator;
    b[4] = exp(-dt / long_s);
    b[5] = exp(-dt / variable_s);
}

static double ms_decay_step(ms_decay_state *state, const double b[6], double input)
{
    double output;
    double memory;
    if (input < state->output) {
        if (state->output > state->memory) {
            memory = state->output * b[0] - state->memory * b[1];
            output = fmax(input, state->output * b[2] - state->memory * b[3]);
            memory = fmin(memory, output);
        } else {
            output = fmax(input, state->output * b[4]);
            memory = output;
        }
    } else {
        output = input;
        memory = state->memory < input
            ? (state->memory - input) * b[5] + input : input;
    }
    state->output = output;
    state->memory = memory;
    return output;
}

ms_status ms_loudness_time_frame_count(
    size_t sample_count, uint32_t sample_rate_hz, size_t *frame_count)
{
    if (frame_count == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (sample_rate_hz != MS_SAMPLE_RATE_HZ) {
        return MS_ERROR_INVALID_SAMPLE_RATE;
    }
    if (sample_count > (size_t)PTRDIFF_MAX / sizeof(double)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    if (sample_count == 0u) {
        return MS_ERROR_INPUT_TOO_SHORT;
    }
    *frame_count = 1u + (sample_count - 1u) / 96u;
    return MS_OK;
}

ms_status ms_loudness_time(
    const double *pressure_pa, size_t sample_count, uint32_t sample_rate_hz,
    ms_sound_field field, double *loudness_sone, size_t loudness_capacity,
    double *specific_sone_per_bark, size_t specific_capacity,
    size_t *frame_count_written)
{
    ms_filter_state filters[MS_THIRD_OCTAVE_BANDS] = {{{{0.0}}}};
    ms_decay_state decay[21] = {{0.0, 0.0}};
    double energy[MS_THIRD_OCTAVE_BANDS][3] = {{0.0}};
    double smoothing[MS_THIRD_OCTAVE_BANDS];
    double previous_core[21] = {0.0};
    double b[6];
    double previous_total = 0.0;
    double fast = 0.0;
    double slow = 0.0;
    const double fast_a = exp(-1.0 / (48000.0 * 0.0035));
    const double slow_a = exp(-1.0 / (48000.0 * 0.070));
    double *temporary;
    double *specific;
    size_t frames, specific_count, total_count, bytes, sample, band;
    ms_status status;
    if (frame_count_written == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    *frame_count_written = 0u;
    status = ms_loudness_time_frame_count(sample_count, sample_rate_hz, &frames);
    if (status != MS_OK) {
        return status;
    }
    if (!ms_is_valid_field(field)) {
        return MS_ERROR_INVALID_SOUND_FIELD;
    }
    if (loudness_sone == NULL || specific_sone_per_bark == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (!ms_checked_multiply(frames, MS_BARK_BANDS, &specific_count)
        || !ms_checked_multiply(frames, MS_BARK_BANDS + 1u, &total_count)
        || !ms_checked_multiply(total_count, sizeof(double), &bytes)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    if (loudness_capacity < frames || specific_capacity < specific_count) {
        *frame_count_written = frames;
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    status = ms_validate_signal(pressure_pa, sample_count, sample_rate_hz);
    if (status != MS_OK) {
        return status;
    }
    temporary = (double *)malloc(bytes);
    if (temporary == NULL) {
        return MS_ERROR_ALLOCATION;
    }
    specific = temporary + frames;
    ms_decay_coefficients(b);
    for (band = 0u; band < MS_THIRD_OCTAVE_BANDS; ++band) {
        const double centre = 1000.0 * pow(10.0, ((double)band - 16.0) / 10.0);
        const double tau = 2.0 / (3.0 * fmin(centre, 1000.0));
        smoothing[band] = exp(-1.0 / (48000.0 * tau));
    }
    for (sample = 0u; sample < sample_count; ++sample) {
        for (band = 0u; band < MS_THIRD_OCTAVE_BANDS; ++band) {
            const double filtered = ms_filter_sample(&filters[band], band, pressure_pa[sample]);
            const double a = smoothing[band];
            double value = filtered * filtered;
            size_t stage;
            if (!ms_is_finite(value)) {
                status = MS_ERROR_NUMERICAL;
                goto cleanup;
            }
            for (stage = 0u; stage < 3u; ++stage) {
                value = (1.0 - a) * value + a * energy[band][stage];
                energy[band][stage] = value;
            }
        }
        if (sample % 24u == 0u) {
            double levels[MS_THIRD_OCTAVE_BANDS];
            double core[21], decayed[21], pattern[MS_BARK_BANDS];
            double total;
            size_t step;
            for (band = 0u; band < MS_THIRD_OCTAVE_BANDS; ++band) {
                levels[band] = ms_energy_level(energy[band][2]);
            }
            ms_core_loudness(levels, field, core);
            for (band = 0u; band < 21u; ++band) {
                if (!ms_is_finite(core[band])) {
                    status = MS_ERROR_NUMERICAL;
                    goto cleanup;
                }
                if (sample != 0u) {
                    for (step = 1u; step < 24u; ++step) {
                        const double input = previous_core[band]
                            + (core[band] - previous_core[band]) * ((double)step / 24.0);
                        (void)ms_decay_step(&decay[band], b, input);
                    }
                }
                decayed[band] = ms_decay_step(&decay[band], b, core[band]);
                previous_core[band] = core[band];
            }
            total = ms_specific_loudness(decayed, pattern);
            if (sample != 0u) {
                for (step = 1u; step < 24u; ++step) {
                    const double input = previous_total
                        + (total - previous_total) * ((double)step / 24.0);
                    fast = (1.0 - fast_a) * input + fast_a * fast;
                    slow = (1.0 - slow_a) * input + slow_a * slow;
                }
            }
            fast = (1.0 - fast_a) * total + fast_a * fast;
            slow = (1.0 - slow_a) * total + slow_a * slow;
            previous_total = total;
            if (!ms_is_finite(total) || !ms_is_finite(fast) || !ms_is_finite(slow)) {
                status = MS_ERROR_NUMERICAL;
                goto cleanup;
            }
            if (sample % 96u == 0u) {
                const size_t frame = sample / 96u;
                temporary[frame] = 0.47 * fast + 0.53 * slow;
                memcpy(specific + frame * MS_BARK_BANDS, pattern, sizeof(pattern));
            }
        }
    }
    memcpy(loudness_sone, temporary, frames * sizeof(double));
    memcpy(specific_sone_per_bark, specific, specific_count * sizeof(double));
    *frame_count_written = frames;
cleanup:
    free(temporary);
    return status;
}
