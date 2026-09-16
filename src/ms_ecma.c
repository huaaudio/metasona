// SPDX-License-Identifier: GPL-3.0-only
// MetaSona author: Jiahua Zhang, 2026.
// C adaptation of SQAT/RefMap ECMA-418-2 Section 5 by Mike JB Lotinga and
// Matt Torjussen, University of Salford (GPL-3.0-only), and MoSQITo by
// Eomys/contributors (Apache-2.0). See NOTICE. Changes: checked C11 buffers,
// per-call filter state, reusable one-band processing, and explicit editions.

#include "ms_ecma.h"

#include <math.h>
#include <stdlib.h>

/* Columns b0,b1,b2,a1,a2. SQAT shmOutMidEarFilter, 2025 coefficients. */
static const double ear2025[8][5] = {
    {1.015896020255593, -1.925298877776079, 0.922118060364679, -1.925298877776079, 0.938014080620272},
    {0.958943219304445, -1.806088011849494, 0.876438777856084, -1.806088011849494, 0.835381997160530},
    {0.961371976333197, -1.763632154338248, 0.821787991845146, -1.763632154338248, 0.783159968178343},
    {2.225803503609735, -1.434650484792157, -0.498204282194628, -1.434650484792157, 0.727599221415107},
    {0.471735128494163, -0.366091796830044, 0.244144703885020, -0.366091796830044, -0.284120167620817},
    {0.115267139824401, 0.0, -0.115267139824401, -1.796002566692014, 0.805837815618546},
    {0.988029297230954, -1.91243380293387, 0.926131550180785, -1.912433802933871, 0.914160847411739},
    {1.952237687301361, 0.162319983017519, -0.667994113035186, 0.162319983017519, 0.284243574266175}
};
static const double ear2022[8][5] = {
    {1.015896, -1.925299, 0.922118, -1.925299, 0.938014},
    {0.958943, -1.806088, 0.876439, -1.806088, 0.835382},
    {0.961372, -1.763632, 0.821788, -1.763632, 0.783160},
    {2.225804, -1.434650, -0.498204, -1.434650, 0.727599},
    {0.471735, -0.366092, 0.244145, -0.366092, -0.284120},
    {0.115267, 0.0, -0.115267, -1.796003, 0.805838},
    {0.988029, -1.912434, 0.926132, -1.912434, 0.914161},
    {1.952238, 0.162320, -0.667994, 0.162320, 0.284244}
};
static const double threshold[MS_ECMA_BANDS] = {
    .3310,.1625,.1051,.0757,.0576,.0453,.0365,.0298,.0247,.0207,
    .0176,.0151,.0131,.0115,.0103,.0093,.0086,.0081,.0077,.0074,
    .0073,.0072,.0071,.0072,.0073,.0074,.0076,.0079,.0082,.0086,
    .0092,.0100,.0109,.0122,.0138,.0157,.0172,.0180,.0180,.0177,
    .0176,.0177,.0182,.0190,.0202,.0217,.0237,.0263,.0296,.0339,
    .0398,.0485,.0622
};

double ms_ecma_centre(size_t band)
{
    return (81.9289 / 0.1618) * sinh(0.1618 * (0.5 * (double)(band + 1u)));
}

double ms_ecma_basis(double rms_pa, size_t band, unsigned edition)
{
    static const double v[9] = {1.0,.6602,.0864,.6384,.0328,.4068,.2082,.3994,.6434};
    double value = .0211668 * rms_pa / MS_REFERENCE_PRESSURE_PA;
    size_t j;
    if (band >= MS_ECMA_BANDS || rms_pa < 0.0 || !ms_is_finite(rms_pa)) return NAN;
    if (edition != 2022u && edition != 2025u) return NAN;
    for (j = 1u; j < 9u; ++j) {
        double p = MS_REFERENCE_PRESSURE_PA * pow(10.0, (5.0 + 10.0 * (double)j) / 20.0);
        value *= pow(1.0 + pow(rms_pa / p, 1.5), (v[j] - v[j - 1u]) / 1.5);
    }
    if (edition == 2025u) value *= 1.00132;
    if (!ms_is_finite(value)) return NAN;
    value -= threshold[band];
    return value > 0.0 ? value : 0.0;
}

ms_status ms_ecma_prepare(const double *input, size_t count, size_t block,
    size_t hop, int pad_end, unsigned edition, ms_sound_field field,
    double **output, size_t *output_count)
{
    size_t padded, rounded, j, stage;
    double *signal;
    ms_status status;
    const double (*coeff)[5];
    if (output == NULL || output_count == NULL) return MS_ERROR_NULL_POINTER;
    *output = NULL;
    *output_count = 0u;
    if (!ms_is_valid_field(field) || block < 2u || hop == 0u || hop > block ||
        (edition != 2022u && edition != 2025u)) return MS_ERROR_INVALID_ARGUMENT;
    if (count < 240u) return MS_ERROR_INPUT_TOO_SHORT;
    if (block > SIZE_MAX - count) return MS_ERROR_SIZE_OVERFLOW;
    padded = count + block;
    if (pad_end) {
        if (padded > SIZE_MAX - (hop - 1u)) return MS_ERROR_SIZE_OVERFLOW;
        rounded = ((padded + hop - 1u) / hop) * hop;
        if (rounded > SIZE_MAX - block) return MS_ERROR_SIZE_OVERFLOW;
        padded = rounded + block;
    }
    if (padded > (size_t)PTRDIFF_MAX / sizeof(double)) return MS_ERROR_SIZE_OVERFLOW;
    status = ms_validate_signal(input, count, MS_SAMPLE_RATE_HZ);
    if (status != MS_OK) return status;
    signal = (double *)calloc(padded, sizeof(double));
    if (signal == NULL) return MS_ERROR_ALLOCATION;
    for (j = 0u; j < count; ++j) {
        double fade = j < 240u ? .5 - .5 * cos(MS_PI * (double)j / 240.0) : 1.0;
        signal[block + j] = input[j] * fade;
    }
    coeff = edition == 2025u ? ear2025 : ear2022;
    for (stage = field == MS_SOUND_FIELD_DIFFUSE ? 2u : 0u; stage < 8u; ++stage) {
        double z1 = 0.0, z2 = 0.0;
        const double *c = coeff[stage];
        for (j = block; j < padded; ++j) {
            double x = signal[j];
            double y = c[0] * x + z1;
            z1 = c[1] * x - c[3] * y + z2;
            z2 = c[2] * x - c[4] * y;
            signal[j] = y;
            if (!ms_is_finite(y)) { free(signal); return MS_ERROR_NUMERICAL; }
        }
    }
    *output = signal;
    *output_count = padded;
    return MS_OK;
}

ms_status ms_ecma_band(const double *input, size_t count, size_t band, double *output)
{
    static const double choose[6] = {1.0,5.0,10.0,10.0,5.0,1.0};
    static const double euler[5] = {0.0,1.0,11.0,11.0,1.0};
    double ar[6], ai[6], br[6], bi[6];
    double zr[5] = {0}, zi[5] = {0};
    double frequency, delta, d, gain, angle;
    size_t j, m;
    if (input == NULL || output == NULL) return MS_ERROR_NULL_POINTER;
    if (band >= MS_ECMA_BANDS) return MS_ERROR_INVALID_ARGUMENT;
    frequency = ms_ecma_centre(band);
    delta = hypot(81.9289, .1618 * frequency);
    d = exp(-delta * 512.0 / (48000.0 * 70.0));
    gain = pow(1.0 - d, 5.0) / (d + 11.0*d*d + 11.0*d*d*d + d*d*d*d);
    angle = 2.0 * MS_PI * frequency / 48000.0;
    for (m = 0u; m < 6u; ++m) {
        double c = cos(angle * (double)m), s = sin(angle * (double)m);
        double a = choose[m] * pow(-d, (double)m);
        double b = m < 5u ? gain * pow(d, (double)m) * euler[m] : 0.0;
        ar[m] = a*c; ai[m] = a*s; br[m] = b*c; bi[m] = b*s;
    }
    /* Complex transposed direct form II, matching the reference filter
     * operation order. No C complex types: portable to MSVC C11. */
    for (j = 0u; j < count; ++j) {
        double x = input[j];
        double yr = zr[0], yi = zi[0]; /* b0 = 0 */
        for (m = 0u; m < 4u; ++m) {
            zr[m] = br[m+1u]*x - (ar[m+1u]*yr - ai[m+1u]*yi) + zr[m+1u];
            zi[m] = bi[m+1u]*x - (ar[m+1u]*yi + ai[m+1u]*yr) + zi[m+1u];
        }
        zr[4] = -(ar[5]*yr - ai[5]*yi);
        zi[4] = -(ar[5]*yi + ai[5]*yr);
        output[j] = 2.0 * yr;
        if (!ms_is_finite(output[j])) return MS_ERROR_NUMERICAL;
    }
    return MS_OK;
}

ms_status ms_ecma_basis_frames(const double *bandpass, size_t count,
    size_t start, size_t block, size_t hop, size_t frames, size_t band,
    unsigned edition, int mosqito_indexing, double *output)
{
    size_t frame, j, last, extent;
    if (bandpass == NULL || output == NULL) return MS_ERROR_NULL_POINTER;
    if (block < 2u || hop == 0u || frames == 0u || band >= MS_ECMA_BANDS ||
        (edition != 2022u && edition != 2025u)) return MS_ERROR_INVALID_ARGUMENT;
    if (!ms_checked_multiply(frames - 1u, hop, &last) || last > SIZE_MAX - start)
        return MS_ERROR_SIZE_OVERFLOW;
    last += start;
    extent = block - (mosqito_indexing ? 0u : 1u);
    if (last >= count || extent >= count - last) return MS_ERROR_INPUT_TOO_SHORT;
    for (frame = 0u; frame < frames; ++frame) {
        double energy = 0.0;
        size_t offset = start + frame * hop;
        for (j = 0u; j < block; ++j) {
            size_t index = offset + j + ((mosqito_indexing && j == block - 1u) ? 1u : 0u);
            double x = bandpass[index];
            if (x > 0.0) energy += x*x;
        }
        if (!ms_is_finite(energy)) return MS_ERROR_NUMERICAL;
        output[frame] = ms_ecma_basis(sqrt(2.0 * energy / (double)block), band, edition);
        if (!ms_is_finite(output[frame])) return MS_ERROR_NUMERICAL;
    }
    return MS_OK;
}
