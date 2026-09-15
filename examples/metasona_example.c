/* SPDX-FileCopyrightText: 2026 MetaSona contributors */
/* SPDX-License-Identifier: Apache-2.0 */

#include <metasona/metasona.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define EXAMPLE_PI 3.14159265358979323846264338327950288

static int report_error(const char *operation, ms_status status)
{
    (void)fprintf(
        stderr,
        "%s failed (%d): %s\n",
        operation,
        (int)status,
        ms_status_string(status));
    return EXIT_FAILURE;
}

int main(void)
{
    const size_t sample_count = MS_SAMPLE_RATE_HZ;
    const double rms_pressure_pa = 0.00002 * 1000.0; /* 60 dB SPL. */
    const double peak_pressure_pa = sqrt(2.0) * rms_pressure_pa;
    double *pressure_pa = NULL;
    double specific[MS_BARK_BANDS];
    double loudness_sone = 0.0;
    double sharpness_acum = 0.0;
    double *time_loudness = NULL;
    double *time_specific = NULL;
    size_t frame_count = 0u;
    size_t frame_count_written = 0u;
    size_t index;
    ms_status status;

    pressure_pa = (double *)malloc(sample_count * sizeof(*pressure_pa));
    if (pressure_pa == NULL) {
        (void)fprintf(stderr, "pressure allocation failed\n");
        return EXIT_FAILURE;
    }
    for (index = 0u; index < sample_count; ++index) {
        const double time_s = (double)index / (double)MS_SAMPLE_RATE_HZ;
        pressure_pa[index] = peak_pressure_pa
                           * sin(2.0 * EXAMPLE_PI * 1000.0 * time_s);
    }

    status = ms_loudness_stationary(
        pressure_pa,
        sample_count,
        MS_SAMPLE_RATE_HZ,
        MS_SOUND_FIELD_FREE,
        &loudness_sone,
        specific,
        MS_BARK_BANDS);
    if (status != MS_OK) {
        free(pressure_pa);
        return report_error("ms_loudness_stationary", status);
    }

    status = ms_sharpness_din(
        specific,
        MS_BARK_BANDS,
        &sharpness_acum);
    if (status != MS_OK) {
        free(pressure_pa);
        return report_error("ms_sharpness_din", status);
    }

    status = ms_loudness_time_frame_count(
        sample_count,
        MS_SAMPLE_RATE_HZ,
        &frame_count);
    if (status != MS_OK) {
        free(pressure_pa);
        return report_error("ms_loudness_time_frame_count", status);
    }
    if (frame_count > SIZE_MAX / MS_BARK_BANDS) {
        free(pressure_pa);
        (void)fprintf(stderr, "specific-output size overflow\n");
        return EXIT_FAILURE;
    }
    time_loudness = (double *)malloc(frame_count * sizeof(*time_loudness));
    time_specific = (double *)malloc(
        frame_count * MS_BARK_BANDS * sizeof(*time_specific));
    if (time_loudness == NULL || time_specific == NULL) {
        free(time_specific);
        free(time_loudness);
        free(pressure_pa);
        (void)fprintf(stderr, "time-output allocation failed\n");
        return EXIT_FAILURE;
    }

    status = ms_loudness_time(
        pressure_pa,
        sample_count,
        MS_SAMPLE_RATE_HZ,
        MS_SOUND_FIELD_FREE,
        time_loudness,
        frame_count,
        time_specific,
        frame_count * MS_BARK_BANDS,
        &frame_count_written);
    if (status != MS_OK) {
        free(time_specific);
        free(time_loudness);
        free(pressure_pa);
        return report_error("ms_loudness_time", status);
    }

    (void)printf("native version: %s (ABI %u)\n",
                 ms_version_string(), (unsigned int)ms_abi_version());
    (void)printf("stationary loudness: %.6f sone\n", loudness_sone);
    (void)printf("sharpness: %.6f acum\n", sharpness_acum);
    (void)printf("time-varying frames: %llu (first %.6f sone)\n",
                 (unsigned long long)frame_count_written, time_loudness[0]);

    free(time_specific);
    free(time_loudness);
    free(pressure_pa);
    return EXIT_SUCCESS;
}
