// SPDX-License-Identifier: Apache-2.0

#include "ms_internal.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MS_AURES_FRAME_SAMPLES 3840u
#define MS_AURES_HOP_SAMPLES 1920u
#define MS_AURES_SPECTRUM_BINS (MS_AURES_FRAME_SAMPLES / 2u + 1u)
#define MS_AURES_MIN_TONE_BIN 2u
#define MS_AURES_MAX_TONE_BIN 400u
#define MS_AURES_MAX_TONES (MS_AURES_MAX_TONE_BIN + 1u)
#define MS_AURES_TONAL_RADIUS 3u
#define MS_AURES_DB_SENTINEL (-10000.0)
#define MS_AURES_LEVEL_FLOOR_DB (-120.0)
#define MS_HANN_HALF_POWER_WIDTH_BINS 1.44092
#define MS_AURES_OFFSET_REFINEMENT_STEPS 32u

typedef struct ms_aures_tone {
    size_t bin;
    double frequency_hz;
    double bark;
    double level_db;
    double intrinsic_width_bark;
    double excess_db;
} ms_aures_tone;

static double ms_add_levels_db(double left_db, double right_db)
{
    double maximum;
    double difference;
    if (left_db <= MS_AURES_DB_SENTINEL) {
        return right_db;
    }
    if (right_db <= MS_AURES_DB_SENTINEL) {
        return left_db;
    }
    maximum = left_db > right_db ? left_db : right_db;
    difference = fabs(left_db - right_db);
    if (difference > 300.0) {
        return maximum;
    }
    return maximum + 10.0 * log10(1.0 + pow(10.0, -0.1 * difference));
}

static double ms_hearing_threshold_db(double frequency_hz)
{
    const double frequency_khz = frequency_hz / 1000.0;
    return 3.64 * pow(frequency_khz, -0.8)
         - 6.5 * exp(-0.6 * (frequency_khz - 3.3)
                          * (frequency_khz - 3.3))
         + 0.001 * pow(frequency_khz, 4.0);
}

static double ms_interpolate_crossing(
    double first_level,
    double second_level,
    double target_level)
{
    const double difference = second_level - first_level;
    double fraction;
    if (fabs(difference) <= DBL_EPSILON) {
        return 0.5;
    }
    fraction = (target_level - first_level) / difference;
    if (fraction < 0.0) {
        return 0.0;
    }
    if (fraction > 1.0) {
        return 1.0;
    }
    return fraction;
}

/*
 * Magnitude of the exact symmetric-Hann transform at an offset expressed in
 * DFT bins.  The three Dirichlet terms follow directly from expanding
 * 0.5 - 0.5*cos(2*pi*n/(N-1)); no sampled correction table is used.
 */
static double ms_hann_response_magnitude(double offset_bins)
{
    const double angular_offset = 2.0 * MS_PI * offset_bins
                                / (double)MS_AURES_FRAME_SAMPLES;
    const double angular_hann = 2.0 * MS_PI
                              / (double)(MS_AURES_FRAME_SAMPLES - 1u);
    const double offsets[3] = {
        angular_offset,
        angular_offset - angular_hann,
        angular_offset + angular_hann
    };
    const double coefficients[3] = {0.5, 0.25, 0.25};
    double response = 0.0;
    size_t index;

    for (index = 0u; index < 3u; ++index) {
        const double denominator = sin(0.5 * offsets[index]);
        double dirichlet;
        if (fabs(denominator) <= DBL_EPSILON) {
            dirichlet = (double)MS_AURES_FRAME_SAMPLES;
        } else {
            dirichlet = sin(
                0.5 * (double)MS_AURES_FRAME_SAMPLES * offsets[index])
                / denominator;
        }
        response += coefficients[index] * dirichlet;
    }
    return fabs(response);
}

static double ms_hann_response_level_db(double offset_bins)
{
    const double coherent_gain =
        0.5 * (double)(MS_AURES_FRAME_SAMPLES - 1u);
    const double magnitude = ms_hann_response_magnitude(offset_bins);
    if (magnitude <= DBL_MIN) {
        return MS_AURES_DB_SENTINEL;
    }
    return 20.0 * log10(magnitude / coherent_gain);
}

static double ms_hann_parabolic_offset(
    double true_offset,
    double *interpolated_peak_db)
{
    const double left = ms_hann_response_level_db(-1.0 - true_offset);
    const double centre = ms_hann_response_level_db(-true_offset);
    const double right = ms_hann_response_level_db(1.0 - true_offset);
    const double denominator = left - 2.0 * centre + right;
    double offset = 0.0;
    double peak = centre;

    if (fabs(denominator) > DBL_EPSILON) {
        offset = 0.5 * (left - right) / denominator;
        if (offset < -0.5) {
            offset = -0.5;
        } else if (offset > 0.5) {
            offset = 0.5;
        }
        peak -= 0.25 * (left - right) * offset;
    }
    if (interpolated_peak_db != NULL) {
        *interpolated_peak_db = peak;
    }
    return offset;
}

/* Undo the small, deterministic offset bias of three-bin parabolic fitting. */
static double ms_hann_true_offset(double observed_offset)
{
    const double sign = observed_offset < 0.0 ? -1.0 : 1.0;
    const double target = fabs(observed_offset);
    double lower = 0.0;
    double upper = 0.5;
    size_t iteration;

    if (target <= DBL_EPSILON) {
        return 0.0;
    }
    if (target >= 0.5) {
        return sign * 0.5;
    }
    for (iteration = 0u;
         iteration < MS_AURES_OFFSET_REFINEMENT_STEPS;
         ++iteration) {
        const double midpoint = 0.5 * (lower + upper);
        const double estimated = fabs(
            ms_hann_parabolic_offset(sign * midpoint, NULL));
        if (estimated < target) {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }
    return sign * 0.5 * (lower + upper);
}

static double ms_tone_intrinsic_width_bark(
    const double bin_levels_db[MS_AURES_SPECTRUM_BINS],
    size_t peak_bin,
    double peak_frequency_hz,
    double peak_level_db)
{
    const double bin_width_hz = (double)MS_SAMPLE_RATE_HZ
                              / (double)MS_AURES_FRAME_SAMPLES;
    const double target_level = peak_level_db - 3.010299956639812;
    size_t lower = peak_bin;
    size_t upper = peak_bin;
    double lower_frequency;
    double upper_frequency;
    double measured_width;
    double window_width;

    while (lower > 0u && bin_levels_db[lower] >= target_level) {
        --lower;
    }
    if (lower == 0u && bin_levels_db[lower] >= target_level) {
        lower_frequency = 0.0;
    } else {
        const double fraction = ms_interpolate_crossing(
            bin_levels_db[lower], bin_levels_db[lower + 1u], target_level);
        lower_frequency = ((double)lower + fraction) * bin_width_hz;
    }

    while (upper + 1u < MS_AURES_SPECTRUM_BINS
           && bin_levels_db[upper] >= target_level) {
        ++upper;
    }
    if (upper + 1u == MS_AURES_SPECTRUM_BINS
        && bin_levels_db[upper] >= target_level) {
        upper_frequency = 0.5 * (double)MS_SAMPLE_RATE_HZ;
    } else {
        const double fraction = ms_interpolate_crossing(
            bin_levels_db[upper - 1u], bin_levels_db[upper], target_level);
        upper_frequency = ((double)(upper - 1u) + fraction) * bin_width_hz;
    }

    measured_width = ms_hz_to_bark(upper_frequency)
                   - ms_hz_to_bark(lower_frequency);
    /*
     * Remove the Hann analysis window's full -3.0103 dB main-lobe width,
     * expressed in DFT-bin units.  Subtracting one bin made the estimated
     * intrinsic bandwidth depend strongly on fractional-bin tone position.
     */
    lower_frequency = peak_frequency_hz
                    - 0.5 * MS_HANN_HALF_POWER_WIDTH_BINS * bin_width_hz;
    if (lower_frequency < 0.0) {
        lower_frequency = 0.0;
    }
    upper_frequency = peak_frequency_hz
                    + 0.5 * MS_HANN_HALF_POWER_WIDTH_BINS * bin_width_hz;
    window_width = ms_hz_to_bark(upper_frequency)
                 - ms_hz_to_bark(lower_frequency);
    if (measured_width <= window_width) {
        return 0.0;
    }
    return measured_width - window_width;
}

static ms_status ms_aures_spectrum(
    const double frame[MS_AURES_FRAME_SAMPLES],
    ms_complex spectrum[MS_AURES_FRAME_SAMPLES],
    double bin_levels_db[MS_AURES_SPECTRUM_BINS],
    double bin_energy_db[MS_AURES_SPECTRUM_BINS],
    int *is_silent)
{
    double maximum_absolute = 0.0;
    double scaled_mean = 0.0;
    double window_sum = 0.0;
    double window_square_sum = 0.0;
    double scale_level_db;
    size_t index;
    ms_status status;

    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        const double absolute = fabs(frame[index]);
        if (absolute > maximum_absolute) {
            maximum_absolute = absolute;
        }
    }
    if (maximum_absolute == 0.0) {
        for (index = 0u; index < MS_AURES_SPECTRUM_BINS; ++index) {
            bin_levels_db[index] = MS_AURES_DB_SENTINEL;
            bin_energy_db[index] = MS_AURES_DB_SENTINEL;
        }
        *is_silent = 1;
        return MS_OK;
    }
    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        scaled_mean += frame[index] / maximum_absolute;
    }
    scaled_mean /= (double)MS_AURES_FRAME_SAMPLES;

    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        const double phase = 2.0 * MS_PI * (double)index
                           / (double)(MS_AURES_FRAME_SAMPLES - 1u);
        const double window = 0.5 - 0.5 * cos(phase);
        spectrum[index].re = (frame[index] / maximum_absolute - scaled_mean)
                           * window;
        spectrum[index].im = 0.0;
        window_sum += window;
        window_square_sum += window * window;
    }
    status = ms_fft_forward(spectrum, MS_AURES_FRAME_SAMPLES);
    if (status != MS_OK) {
        return status;
    }

    scale_level_db = 20.0 * log10(maximum_absolute)
                   - 20.0 * log10(MS_REFERENCE_PRESSURE_PA);
    for (index = 0u; index < MS_AURES_SPECTRUM_BINS; ++index) {
        const double magnitude = hypot(spectrum[index].re, spectrum[index].im);
        const double one_sided_factor = (index == 0u
            || index + 1u == MS_AURES_SPECTRUM_BINS) ? 1.0 : 2.0;
        if (!ms_is_finite(magnitude)) {
            return MS_ERROR_NUMERICAL;
        }
        if (magnitude == 0.0) {
            bin_levels_db[index] = MS_AURES_DB_SENTINEL;
            bin_energy_db[index] = MS_AURES_DB_SENTINEL;
        } else {
            bin_levels_db[index] = 20.0 * log10(magnitude)
                + 20.0 * log10(sqrt(one_sided_factor) / window_sum)
                + scale_level_db;
            bin_energy_db[index] = 20.0 * log10(magnitude)
                + 10.0 * log10(one_sided_factor)
                - 10.0 * log10((double)MS_AURES_FRAME_SAMPLES
                              * window_square_sum)
                + scale_level_db;
            if (!ms_is_finite(bin_levels_db[index])
                || !ms_is_finite(bin_energy_db[index])) {
                return MS_ERROR_NUMERICAL;
            }
        }
    }
    *is_silent = 0;
    return MS_OK;
}

static size_t ms_find_aures_tones(
    const double bin_levels_db[MS_AURES_SPECTRUM_BINS],
    ms_aures_tone tones[MS_AURES_MAX_TONES])
{
    const double bin_width_hz = (double)MS_SAMPLE_RATE_HZ
                              / (double)MS_AURES_FRAME_SAMPLES;
    size_t tone_count = 0u;
    size_t bin;
    for (bin = MS_AURES_MIN_TONE_BIN + 3u;
         bin + 3u <= MS_AURES_MAX_TONE_BIN;
         ++bin) {
        const double level = bin_levels_db[bin];
        if (level <= bin_levels_db[bin - 1u]
            || level < bin_levels_db[bin + 1u]
            || level - bin_levels_db[bin - 2u] < 7.0
            || level - bin_levels_db[bin - 3u] < 7.0
            || level - bin_levels_db[bin + 2u] < 7.0
            || level - bin_levels_db[bin + 3u] < 7.0) {
            continue;
        }
        if (tone_count < MS_AURES_MAX_TONES) {
            const double denominator = bin_levels_db[bin - 1u]
                                     - 2.0 * level
                                     + bin_levels_db[bin + 1u];
            double observed_offset = 0.0;
            double corrected_offset;
            double parabolic_bias_db;
            double interpolated_level = level;
            if (fabs(denominator) > DBL_EPSILON) {
                observed_offset = 0.5 * (bin_levels_db[bin - 1u]
                                       - bin_levels_db[bin + 1u])
                                  / denominator;
                if (observed_offset < -0.5) {
                    observed_offset = -0.5;
                } else if (observed_offset > 0.5) {
                    observed_offset = 0.5;
                }
                interpolated_level -= 0.25
                    * (bin_levels_db[bin - 1u] - bin_levels_db[bin + 1u])
                    * observed_offset;
            }
            corrected_offset = ms_hann_true_offset(observed_offset);
            (void)ms_hann_parabolic_offset(
                corrected_offset, &parabolic_bias_db);
            tones[tone_count].bin = bin;
            tones[tone_count].frequency_hz = ((double)bin + corrected_offset)
                                            * bin_width_hz;
            tones[tone_count].bark = ms_hz_to_bark(
                tones[tone_count].frequency_hz);
            tones[tone_count].level_db = interpolated_level
                                       - parabolic_bias_db;
            tones[tone_count].intrinsic_width_bark =
                ms_tone_intrinsic_width_bark(
                    bin_levels_db,
                    bin,
                    tones[tone_count].frequency_hz,
                    interpolated_level);
            tones[tone_count].excess_db = 0.0;
            ++tone_count;
        }
    }
    return tone_count;
}

static void ms_mark_tonal_bins(
    const ms_aures_tone *tones,
    size_t tone_count,
    unsigned char tonal_bins[MS_AURES_SPECTRUM_BINS])
{
    size_t tone_index;
    memset(tonal_bins, 0, MS_AURES_SPECTRUM_BINS * sizeof(*tonal_bins));
    for (tone_index = 0u; tone_index < tone_count; ++tone_index) {
        size_t first;
        size_t last;
        size_t bin;
        if (tones[tone_index].excess_db <= 0.0) {
            continue;
        }
        first = tones[tone_index].bin > MS_AURES_TONAL_RADIUS
            ? tones[tone_index].bin - MS_AURES_TONAL_RADIUS : 0u;
        last = tones[tone_index].bin + MS_AURES_TONAL_RADIUS;
        if (last >= MS_AURES_SPECTRUM_BINS) {
            last = MS_AURES_SPECTRUM_BINS - 1u;
        }
        for (bin = first; bin <= last; ++bin) {
            tonal_bins[bin] = 1u;
        }
    }
}

static double ms_secondary_masking_level_db(
    const ms_aures_tone *tones,
    size_t tone_count,
    size_t selected)
{
    double maximum_level = MS_AURES_DB_SENTINEL;
    double scaled_sum = 0.0;
    size_t other;
    for (other = 0u; other < tone_count; ++other) {
        double slope;
        double masked_level;
        if (other == selected) {
            continue;
        }
        if (tones[other].frequency_hz < tones[selected].frequency_hz) {
            slope = -24.0 - 230.0 / tones[other].frequency_hz
                  + 0.2 * tones[other].level_db;
        } else {
            slope = 27.0;
        }
        masked_level = tones[other].level_db
            - slope * (tones[other].bark - tones[selected].bark);
        if (masked_level > maximum_level) {
            maximum_level = masked_level;
        }
    }
    if (maximum_level <= MS_AURES_DB_SENTINEL) {
        return MS_AURES_DB_SENTINEL;
    }
    for (other = 0u; other < tone_count; ++other) {
        double slope;
        double masked_level;
        if (other == selected) {
            continue;
        }
        if (tones[other].frequency_hz < tones[selected].frequency_hz) {
            slope = -24.0 - 230.0 / tones[other].frequency_hz
                  + 0.2 * tones[other].level_db;
        } else {
            slope = 27.0;
        }
        masked_level = tones[other].level_db
            - slope * (tones[other].bark - tones[selected].bark);
        if (maximum_level - masked_level <= 6000.0) {
            scaled_sum += pow(10.0, (masked_level - maximum_level) / 20.0);
        }
    }
    return maximum_level + 20.0 * log10(scaled_sum);
}

static size_t ms_compute_tone_excesses(
    ms_aures_tone *tones,
    size_t tone_count,
    const double bin_energy_db[MS_AURES_SPECTRUM_BINS])
{
    const double bin_width_hz = (double)MS_SAMPLE_RATE_HZ
                              / (double)MS_AURES_FRAME_SAMPLES;
    unsigned char candidate_bins[MS_AURES_SPECTRUM_BINS];
    size_t audible_count = 0u;
    size_t tone_index;

    memset(candidate_bins, 0, sizeof(candidate_bins));
    for (tone_index = 0u; tone_index < tone_count; ++tone_index) {
        size_t first = tones[tone_index].bin > MS_AURES_TONAL_RADIUS
            ? tones[tone_index].bin - MS_AURES_TONAL_RADIUS : 0u;
        size_t last = tones[tone_index].bin + MS_AURES_TONAL_RADIUS;
        size_t bin;
        if (last >= MS_AURES_SPECTRUM_BINS) {
            last = MS_AURES_SPECTRUM_BINS - 1u;
        }
        for (bin = first; bin <= last; ++bin) {
            candidate_bins[bin] = 1u;
        }
    }

    for (tone_index = 0u; tone_index < tone_count; ++tone_index) {
        double noise_level = MS_AURES_DB_SENTINEL;
        double combined_masking;
        double secondary_level;
        size_t bin;
        for (bin = MS_AURES_MIN_TONE_BIN;
             bin <= MS_AURES_MAX_TONE_BIN;
             ++bin) {
            const double frequency = (double)bin * bin_width_hz;
            const double bark = ms_hz_to_bark(frequency);
            if (candidate_bins[bin] == 0u
                && bark >= tones[tone_index].bark - 0.5
                && bark <= tones[tone_index].bark + 0.5) {
                noise_level = ms_add_levels_db(
                    noise_level, bin_energy_db[bin]);
            }
        }
        combined_masking = ms_add_levels_db(
            noise_level,
            ms_hearing_threshold_db(tones[tone_index].frequency_hz));
        secondary_level = ms_secondary_masking_level_db(
            tones, tone_count, tone_index);
        combined_masking = ms_add_levels_db(combined_masking, secondary_level);
        tones[tone_index].excess_db = tones[tone_index].level_db
                                    - combined_masking;
        if (tones[tone_index].excess_db > 0.0) {
            ++audible_count;
        } else {
            tones[tone_index].excess_db = 0.0;
        }
    }
    return audible_count;
}

static void ms_third_octave_levels_from_spectrum(
    const double bin_energy_db[MS_AURES_SPECTRUM_BINS],
    const unsigned char *excluded_bins,
    double levels_db[MS_THIRD_OCTAVE_BANDS])
{
    const double bin_width_hz = (double)MS_SAMPLE_RATE_HZ
                              / (double)MS_AURES_FRAME_SAMPLES;
    const double sixth_octave = pow(2.0, 1.0 / 6.0);
    const double *centres = ms_third_octave_centres_hz();
    size_t band;
    for (band = 0u; band < MS_THIRD_OCTAVE_BANDS; ++band) {
        const double lower = centres[band] / sixth_octave;
        const double upper = centres[band] * sixth_octave;
        double level = MS_AURES_DB_SENTINEL;
        size_t bin;
        for (bin = 1u; bin < MS_AURES_SPECTRUM_BINS; ++bin) {
            const double frequency = (double)bin * bin_width_hz;
            if (frequency >= lower && frequency < upper
                && (excluded_bins == NULL || excluded_bins[bin] == 0u)) {
                level = ms_add_levels_db(level, bin_energy_db[bin]);
            }
        }
        levels_db[band] = level <= MS_AURES_DB_SENTINEL
            ? MS_AURES_LEVEL_FLOOR_DB : level;
        if (levels_db[band] < MS_AURES_LEVEL_FLOOR_DB) {
            levels_db[band] = MS_AURES_LEVEL_FLOOR_DB;
        }
    }
}

/*
 * Remove fitted narrowband sinusoids before estimating residual-noise
 * loudness.  This compensates Hann leakage outside the seven-line tonal core
 * without enlarging that model-defined core.  The fit has fixed O(N*T)
 * cost, where T is already bounded by MS_AURES_MAX_TONES, and reuses the
 * caller's FFT workspace.
 */
static ms_status ms_aures_residual_energy(
    const double frame[MS_AURES_FRAME_SAMPLES],
    const ms_aures_tone *tones,
    size_t tone_count,
    ms_complex spectrum[MS_AURES_FRAME_SAMPLES],
    double residual_energy_db[MS_AURES_SPECTRUM_BINS])
{
    double maximum_absolute = 0.0;
    double scaled_mean = 0.0;
    double window_square_sum = 0.0;
    double scale_level_db;
    size_t index;
    size_t tone_index;
    ms_status status;

    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        const double absolute = fabs(frame[index]);
        if (absolute > maximum_absolute) {
            maximum_absolute = absolute;
        }
    }
    if (maximum_absolute == 0.0) {
        for (index = 0u; index < MS_AURES_SPECTRUM_BINS; ++index) {
            residual_energy_db[index] = MS_AURES_DB_SENTINEL;
        }
        return MS_OK;
    }
    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        scaled_mean += frame[index] / maximum_absolute;
    }
    scaled_mean /= (double)MS_AURES_FRAME_SAMPLES;

    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        const double phase = 2.0 * MS_PI * (double)index
                           / (double)(MS_AURES_FRAME_SAMPLES - 1u);
        const double window = 0.5 - 0.5 * cos(phase);
        spectrum[index].re = frame[index] / maximum_absolute - scaled_mean;
        spectrum[index].im = window;
        window_square_sum += window * window;
    }

    for (tone_index = 0u; tone_index < tone_count; ++tone_index) {
        const double angular_step = 2.0 * MS_PI
                                  * tones[tone_index].frequency_hz
                                  / (double)MS_SAMPLE_RATE_HZ;
        const double half_step_sine = sin(0.5 * angular_step);
        double mean_cosine = 0.0;
        double mean_sine = 0.0;
        double cosine = 1.0;
        double sine = 0.0;
        const double cosine_step = cos(angular_step);
        const double sine_step = sin(angular_step);
        double cosine_square = 0.0;
        double sine_square = 0.0;
        double cross = 0.0;
        double cosine_signal = 0.0;
        double sine_signal = 0.0;
        double determinant;
        double cosine_amplitude;
        double sine_amplitude;

        if (tones[tone_index].excess_db <= 0.0) {
            continue;
        }
        if (fabs(half_step_sine) > DBL_EPSILON) {
            const double sum_scale = sin(
                0.5 * (double)MS_AURES_FRAME_SAMPLES * angular_step)
                / ((double)MS_AURES_FRAME_SAMPLES * half_step_sine);
            const double mean_phase = 0.5
                * (double)(MS_AURES_FRAME_SAMPLES - 1u) * angular_step;
            mean_cosine = sum_scale * cos(mean_phase);
            mean_sine = sum_scale * sin(mean_phase);
        }

        for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
            const double cosine_basis = cosine - mean_cosine;
            const double sine_basis = sine - mean_sine;
            const double weight = spectrum[index].im * spectrum[index].im;
            const double next_cosine = cosine * cosine_step
                                     - sine * sine_step;
            const double next_sine = sine * cosine_step
                                   + cosine * sine_step;
            cosine_square += weight * cosine_basis * cosine_basis;
            sine_square += weight * sine_basis * sine_basis;
            cross += weight * cosine_basis * sine_basis;
            cosine_signal += weight * cosine_basis * spectrum[index].re;
            sine_signal += weight * sine_basis * spectrum[index].re;
            cosine = next_cosine;
            sine = next_sine;
        }
        determinant = cosine_square * sine_square - cross * cross;
        if (determinant <= DBL_EPSILON * cosine_square * sine_square) {
            continue;
        }
        cosine_amplitude = (cosine_signal * sine_square
                          - sine_signal * cross) / determinant;
        sine_amplitude = (sine_signal * cosine_square
                        - cosine_signal * cross) / determinant;
        if (!ms_is_finite(cosine_amplitude)
            || !ms_is_finite(sine_amplitude)) {
            return MS_ERROR_NUMERICAL;
        }

        cosine = 1.0;
        sine = 0.0;
        for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
            const double next_cosine = cosine * cosine_step
                                     - sine * sine_step;
            const double next_sine = sine * cosine_step
                                   + cosine * sine_step;
            spectrum[index].re -= cosine_amplitude
                * (cosine - mean_cosine)
                + sine_amplitude * (sine - mean_sine);
            cosine = next_cosine;
            sine = next_sine;
        }
    }

    for (index = 0u; index < MS_AURES_FRAME_SAMPLES; ++index) {
        spectrum[index].re *= spectrum[index].im;
        spectrum[index].im = 0.0;
    }
    status = ms_fft_forward(spectrum, MS_AURES_FRAME_SAMPLES);
    if (status != MS_OK) {
        return status;
    }
    scale_level_db = 20.0 * log10(maximum_absolute)
                   - 20.0 * log10(MS_REFERENCE_PRESSURE_PA);
    for (index = 0u; index < MS_AURES_SPECTRUM_BINS; ++index) {
        const double magnitude = hypot(spectrum[index].re, spectrum[index].im);
        const double one_sided_factor = (index == 0u
            || index + 1u == MS_AURES_SPECTRUM_BINS) ? 1.0 : 2.0;
        if (!ms_is_finite(magnitude)) {
            return MS_ERROR_NUMERICAL;
        }
        if (magnitude == 0.0) {
            residual_energy_db[index] = MS_AURES_DB_SENTINEL;
        } else {
            residual_energy_db[index] = 20.0 * log10(magnitude)
                + 10.0 * log10(one_sided_factor)
                - 10.0 * log10((double)MS_AURES_FRAME_SAMPLES
                              * window_square_sum)
                + scale_level_db;
            if (!ms_is_finite(residual_energy_db[index])) {
                return MS_ERROR_NUMERICAL;
            }
        }
    }
    return MS_OK;
}

static ms_status ms_aures_loudness_weight(
    const double frame[MS_AURES_FRAME_SAMPLES],
    ms_complex spectrum[MS_AURES_FRAME_SAMPLES],
    const double bin_energy_db[MS_AURES_SPECTRUM_BINS],
    const ms_aures_tone *tones,
    size_t tone_count,
    ms_sound_field field,
    double *weight)
{
    unsigned char tonal_bins[MS_AURES_SPECTRUM_BINS];
    double residual_energy_db[MS_AURES_SPECTRUM_BINS];
    double total_levels[MS_THIRD_OCTAVE_BANDS];
    double noise_levels[MS_THIRD_OCTAVE_BANDS];
    double total_loudness;
    double noise_loudness;
    ms_status status;

    status = ms_aures_residual_energy(
        frame, tones, tone_count, spectrum, residual_energy_db);
    if (status != MS_OK) {
        return status;
    }
    ms_mark_tonal_bins(tones, tone_count, tonal_bins);
    ms_third_octave_levels_from_spectrum(bin_energy_db, NULL, total_levels);
    ms_third_octave_levels_from_spectrum(
        residual_energy_db, tonal_bins, noise_levels);
    status = ms_loudness_from_levels(
        total_levels,
        MS_THIRD_OCTAVE_BANDS,
        field,
        &total_loudness,
        NULL,
        0u);
    if (status != MS_OK) {
        return status;
    }
    status = ms_loudness_from_levels(
        noise_levels,
        MS_THIRD_OCTAVE_BANDS,
        field,
        &noise_loudness,
        NULL,
        0u);
    if (status != MS_OK) {
        return status;
    }
    if (total_loudness <= DBL_EPSILON) {
        *weight = 0.0;
        return MS_OK;
    }
    *weight = 1.0 - noise_loudness / total_loudness;
    if (*weight < 0.0) {
        *weight = 0.0;
    } else if (*weight > 1.0) {
        *weight = 1.0;
    }
    return MS_OK;
}

static double ms_aures_tonal_weight(
    const ms_aures_tone *tones,
    size_t tone_count)
{
    const double aures_exponent = 0.29;
    double squared_sum = 0.0;
    size_t tone_index;
    for (tone_index = 0u; tone_index < tone_count; ++tone_index) {
        double bandwidth_weight;
        double frequency_weight;
        double level_weight;
        double primed_weight;
        const double frequency = tones[tone_index].frequency_hz;
        if (tones[tone_index].excess_db <= 0.0) {
            continue;
        }
        bandwidth_weight = 0.13
            / (tones[tone_index].intrinsic_width_bark + 0.13);
        frequency_weight = 1.0 / sqrt(
            1.0 + 0.2 * (frequency / 700.0 + 700.0 / frequency)
                      * (frequency / 700.0 + 700.0 / frequency));
        level_weight = 1.0 - exp(-tones[tone_index].excess_db / 15.0);
        primed_weight = pow(bandwidth_weight, 1.0 / aures_exponent)
                       * frequency_weight * level_weight;
        squared_sum += primed_weight * primed_weight;
    }
    if (squared_sum <= 0.0) {
        return 0.0;
    }
    return pow(sqrt(squared_sum), aures_exponent);
}

static ms_status ms_aures_frame(
    const double frame[MS_AURES_FRAME_SAMPLES],
    ms_sound_field field,
    ms_complex spectrum[MS_AURES_FRAME_SAMPLES],
    double *tonality)
{
    double bin_levels_db[MS_AURES_SPECTRUM_BINS];
    double bin_energy_db[MS_AURES_SPECTRUM_BINS];
    ms_aures_tone tones[MS_AURES_MAX_TONES];
    double loudness_weight;
    double tonal_weight;
    size_t tone_count;
    int is_silent;
    ms_status status = ms_aures_spectrum(
        frame, spectrum, bin_levels_db, bin_energy_db, &is_silent);
    if (status != MS_OK) {
        return status;
    }
    if (is_silent != 0) {
        *tonality = 0.0;
        return MS_OK;
    }
    tone_count = ms_find_aures_tones(bin_levels_db, tones);
    if (tone_count == 0u
        || ms_compute_tone_excesses(tones, tone_count, bin_energy_db) == 0u) {
        *tonality = 0.0;
        return MS_OK;
    }
    status = ms_aures_loudness_weight(
        frame,
        spectrum,
        bin_energy_db,
        tones,
        tone_count,
        field,
        &loudness_weight);
    if (status != MS_OK) {
        return status;
    }
    tonal_weight = ms_aures_tonal_weight(tones, tone_count);
    *tonality = 1.09 * tonal_weight * pow(loudness_weight, 0.79);
    if (!ms_is_finite(*tonality) || *tonality < 0.0) {
        return MS_ERROR_NUMERICAL;
    }
    return MS_OK;
}

ms_status ms_tonality_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count)
{
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
    if (sample_count < MS_AURES_FRAME_SAMPLES) {
        return MS_ERROR_INPUT_TOO_SHORT;
    }
    result = 1u + (sample_count - MS_AURES_FRAME_SAMPLES)
                  / MS_AURES_HOP_SAMPLES;
    *frame_count = result;
    return MS_OK;
}

ms_status ms_tonality_aures(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *tonality,
    size_t tonality_capacity,
    size_t *frame_count_written)
{
    ms_complex *spectrum = NULL;
    double *temporary = NULL;
    size_t frame_count;
    size_t spectrum_bytes;
    size_t output_bytes;
    size_t frame;
    ms_status status;

    if (pressure_pa == NULL || tonality == NULL || frame_count_written == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    status = ms_tonality_frame_count(
        sample_count, sample_rate_hz, &frame_count);
    if (status != MS_OK) {
        return status;
    }
    if (!ms_is_valid_field(field)) {
        return MS_ERROR_INVALID_SOUND_FIELD;
    }
    if (tonality_capacity < frame_count) {
        return MS_ERROR_OUTPUT_TOO_SMALL;
    }
    status = ms_validate_signal(pressure_pa, sample_count, sample_rate_hz);
    if (status != MS_OK) {
        return status;
    }
    if (!ms_checked_multiply(
            MS_AURES_FRAME_SAMPLES, sizeof(*spectrum), &spectrum_bytes)
        || !ms_checked_multiply(frame_count, sizeof(*temporary), &output_bytes)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    spectrum = (ms_complex *)calloc(
        MS_AURES_FRAME_SAMPLES, sizeof(*spectrum));
    temporary = (double *)calloc(frame_count, sizeof(*temporary));
    if (spectrum == NULL || temporary == NULL) {
        free(spectrum);
        free(temporary);
        return MS_ERROR_ALLOCATION;
    }

    for (frame = 0u; frame < frame_count; ++frame) {
        const size_t offset = frame * MS_AURES_HOP_SAMPLES;
        status = ms_aures_frame(
            pressure_pa + offset, field, spectrum, &temporary[frame]);
        if (status != MS_OK) {
            free(spectrum);
            free(temporary);
            return status;
        }
    }
    memcpy(tonality, temporary, output_bytes);
    *frame_count_written = frame_count;
    free(spectrum);
    free(temporary);
    return MS_OK;
}

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
