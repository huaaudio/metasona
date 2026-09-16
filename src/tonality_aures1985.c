// SPDX-License-Identifier: GPL-3.0-only AND Apache-2.0 AND BSD-3-Clause
// Based on the supplied C port; aligned to SQAT e6228b789fc9.
// MetaSona adaptation: Jiahua Zhang, September 2026.
// Changes: MetaSona loudness, complete frames, checked FFT/numeric errors.
// Provenance and notices: THIRD_PARTY.md, LICENSES/BSD-3-Clause-PA-Tonality.txt.

/************************************************************************/
/*  Tonality calculation according to Aures (1985)                      */
/*                                                                      */
/*  Ported from archive/PA/Tonality_Aures1985/Tonality_Aures1985.m.     */
/************************************************************************/

#define TONALITY_AURES1985_BUILD
#include "tonality_aures1985.h"

#include "ms_internal.h"
#include "pocketfft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define AURES_TIME_RESOLUTION 0.250
#define AURES_MIN_FREQUENCY   20.0
#define AURES_MAX_FREQUENCY   5000.0
#define AURES_TONE_THRESHOLD  7.0
/* SQAT logarithm floor; it does not change finite audible levels. */
#define AURES_TINY_VALUE      1e-99
#define AURES_I_REF           4e-10
/* SQAT implementation calibration (upstream gives ideal-model C = 1.09).
 * This is inherited normalization, not a universal Aures constant. */
#define AURES_REFERENCE_SCALE 1.1055
/* SQAT half-power bandwidth threshold in dB, separate from notch width. */
#define AURES_HALF_POWER_DB 3.0
/* Reference skips 5% of each frame when averaging loudness-filter energy. */
#define AURES_LOUDNESS_SKIP_FRACTION 0.05

static int round_to_int(double x)
{
    return (int)floor(x + 0.5);
}

static double complex_abs_at(const double *fftData, int idx)
{
    double re = fftData[2 * idx];
    double im = fftData[2 * idx + 1];
    return sqrt(re * re + im * im);
}

static void hann_window(double *window, int n)
{
    if (n <= 1) {
        if (n == 1) window[0] = 1.0;
        return;
    }

    for (int i = 0; i < n; i++) {
        window[i] = 0.5 - 0.5 * cos(2.0 * M_PI * (double)i / (double)(n - 1));
    }
}

static double mean_array(const double *x, int n)
{
    double s = 0.0;
    for (int i = 0; i < n; i++) s += x[i];
    return (n > 0) ? s / (double)n : 0.0;
}

static double fq_to_bark(double f)
{
    double fk = f / 1000.0;
    return 13.0 * atan(0.76 * fk) + 3.5 * atan((fk / 7.5) * (fk / 7.5));
}

static double threshold_hearing(double f)
{
    double fk = f / 1000.0;
    if (fk <= 0.0) return 130.0;
    return 3.64 * pow(fk, -0.8)
         - 6.5 * exp(-0.6 * (fk - 3.3) * (fk - 3.3))
         + 1e-3 * fk * fk * fk * fk;
}

static double interp_xy(double x0, double y0, double x1, double y1, double x)
{
    /* A flat segment cannot define a unique crossing. The caller reports
     * a numerical error instead of inventing a midpoint. */
    if (x1 == x0) return NAN;
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

static int do_fft(double *data, int n, cfft_plan plan)
{
    (void)n;
    return cfft_forward(plan, data, 1.0);
}

static int do_ifft(double *data, int n, cfft_plan plan)
{
    return cfft_backward(plan, data, 1.0 / (double)n);
}

static int compute_stationary_loudness(const double *signal, int n, double fs,
                                       int soundField, double timeSkip,
                                       double *outLoudness)
{
    double levels[MS_THIRD_OCTAVE_BANDS], specific[MS_BARK_BANDS];
    ms_status status = ms_levels_from_signal_skip_48k(
        signal, (size_t)n, (size_t)(timeSkip * fs), levels);
    if (status == MS_OK) status = ms_loudness_from_levels(
        levels, MS_THIRD_OCTAVE_BANDS, (ms_sound_field)soundField,
        outLoudness, specific, MS_BARK_BANDS);
    return status == MS_OK ? 0 : TonalityAures1985ErrorLoudnessFailed;
}

static int compare_double_ascending(const void *a, const void *b)
{
    double da = *(const double*)a;
    double db = *(const double*)b;
    return (da > db) - (da < db);
}

static double percentile_sorted(const double *sorted, int n, double q)
{
    if (n <= 0) return 0.0;
    if (n == 1) return sorted[0];
    if (q <= 0.0) return sorted[0];
    if (q >= 1.0) return sorted[n - 1];

    double pos = q * (double)(n - 1);
    int lo = (int)floor(pos);
    int hi = lo + 1;
    double frac = pos - (double)lo;
    if (hi >= n) return sorted[n - 1];
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

static int fill_stats(const double *values, const double *time, int n,
                      double timeSkip, double *outStats)
{
    int start = 0;
    double best = fabs(time[0] - timeSkip);

    for (int i = 1; i < n; i++) {
        double d = fabs(time[i] - timeSkip);
        if (d < best) {
            best = d;
            start = i;
        }
    }

    int m = n - start;
    if (m <= 0) return TonalityAures1985ErrorSignalTooShort;

    double mean = 0.0;
    double maxv = values[start];
    double minv = values[start];
    for (int i = start; i < n; i++) {
        double v = values[i];
        mean += v;
        if (v > maxv) maxv = v;
        if (v < minv) minv = v;
    }
    mean /= (double)m;

    double ss = 0.0;
    for (int i = start; i < n; i++) {
        double d = values[i] - mean;
        ss += d * d;
    }
    double stdv = (m > 1) ? sqrt(ss / (double)(m - 1)) : 0.0;

    double *tmp = (double*)calloc((size_t)m, sizeof(double));
    if (!tmp) return TonalityAures1985ErrorMemoryAlloc;
    memcpy(tmp, values + start, (size_t)m * sizeof(double));
    qsort(tmp, (size_t)m, sizeof(double), compare_double_ascending);
    double p5 = percentile_sorted(tmp, m, 0.95);
    free(tmp);

    outStats[0] = mean;
    outStats[1] = stdv;
    outStats[2] = maxv;
    outStats[3] = minv;
    outStats[4] = p5;

    return 0;
}

static int estimate_bandwidths(const double *spl, int n, double fs,
                                const int *toneBins, const double *toneLevels,
                                int nTones, double *outBw)
{
    for (int i = 0; i < nTones; i++) {
        int idx = toneBins[i];
        double target = toneLevels[i] - AURES_HALF_POWER_DB;
        int lowIdx = -1;
        int highIdx = -1;

        for (int j = 0; j <= idx; j++) {
            if (spl[j] < target) lowIdx = j;
        }
        if (lowIdx < 3) lowIdx = 3;
        if (lowIdx + 1 >= n) lowIdx = n - 2;

        for (int j = idx + 1; j < n; j++) {
            if (spl[j] < target) {
                highIdx = j;
                break;
            }
        }
        if (highIdx < 1) highIdx = idx + 1;
        if (highIdx >= n) highIdx = n - 1;

        double fLow = interp_xy(spl[lowIdx], (double)lowIdx * fs / (double)n,
                                spl[lowIdx + 1], (double)(lowIdx + 1) * fs / (double)n,
                                target);
        double fHigh = interp_xy(spl[highIdx - 1], (double)(highIdx - 1) * fs / (double)n,
                                 spl[highIdx], (double)highIdx * fs / (double)n,
                                 target);
        double bw = fHigh - fLow;
        /* Do not fabricate a 1 Hz bandwidth on failed interpolation. */
        if (!isfinite(bw) || bw == 0.0) return -1;
        if (bw < 0.0) return -1;
        outBw[i] = bw;
    }
    return 0;
}

static void remove_tones_preserve_phase(double *spectrum, int n, double fs,
                                        const double *toneFreq,
                                        const double *toneBw,
                                        int nTones)
{
    int singleLen = n / 2 + 1;

    for (int i = 0; i < nTones; i++) {
        double fLow = toneFreq[i] - 0.5 * toneBw[i];
        double fHigh = toneFreq[i] + 0.5 * toneBw[i];
        int low = 0;
        int up = singleLen - 1;

        for (int k = 0; k < singleLen; k++) {
            if ((double)k * fs / (double)n >= fLow) {
                low = k;
                break;
            }
        }
        for (int k = 0; k < singleLen; k++) {
            if ((double)k * fs / (double)n >= fHigh) {
                up = k;
                break;
            }
        }

        if (low < 0) low = 0;
        if (up < low) up = low;
        if (up >= singleLen) up = singleLen - 1;

        int left = (low == 0) ? low : low - 1;
        int right = (up + 1 < singleLen) ? up + 1 : up;
        double mag = 0.5 * (complex_abs_at(spectrum, left) +
                            complex_abs_at(spectrum, right));

        for (int k = low; k <= up; k++) {
            double phase = atan2(spectrum[2 * k + 1], spectrum[2 * k]);
            spectrum[2 * k] = mag * cos(phase);
            spectrum[2 * k + 1] = mag * sin(phase);
        }
    }

    for (int k = 1; k < singleLen - 1; k++) {
        int dst = n - k;
        spectrum[2 * dst] = spectrum[2 * k];
        spectrum[2 * dst + 1] = -spectrum[2 * k + 1];
    }
    spectrum[1] = 0.0;
    if ((n % 2) == 0) spectrum[2 * (n / 2) + 1] = 0.0;
}

static int spl_excess(const double *splCrop, int cropLen, int minIdx,
                       double fs, int nFft, const double *toneFreq,
                       const double *toneLevel, const int *toneCropIdx,
                       int nTones, double *outLx)
{
    for (int i = 0; i < nTones; i++) {
        double toneBark = fq_to_bark(toneFreq[i]);
        double lowBark = toneBark - 0.5;
        double highBark = toneBark + 0.5;
        double egr = 0.0;

        for (int k = 0; k < cropLen; k++) {
            double freq = (double)(minIdx + k) * fs / (double)nFft;
            double bark = fq_to_bark(freq);
            int skipTone = toneCropIdx[i] >= 0
                && (k >= toneCropIdx[i] - 2 && k <= toneCropIdx[i] + 2);
            if (!skipTone && bark >= lowBark && bark <= highBark) {
                /* Dimensionless residual-noise intensity, as in SQAT. */
                egr += pow(10.0, splCrop[k] / 10.0);
            }
        }

        double sumlo = AURES_TINY_VALUE;
        double sumhi = AURES_TINY_VALUE;
        for (int j = 0; j < nTones; j++) {
            if (j == i) continue;

            double barkJ = fq_to_bark(toneFreq[j]);
            double lji;
            if (j < i) {
                double s = -24.0 - (230.0 / toneFreq[j]) + (0.2 * toneLevel[j]);
                lji = toneLevel[j] - s * (barkJ - toneBark);
                sumlo += pow(10.0, lji / 20.0);
            } else {
                lji = toneLevel[j] - 27.0 * (barkJ - toneBark);
                sumhi += pow(10.0, lji / 20.0);
            }
        }

        double aek = sumlo + sumhi;
        double ehs = pow(10.0, threshold_hearing(toneFreq[i]) / 10.0);
        double den = (nTones == 1) ? (egr + ehs) : (aek * aek + egr + ehs);
        double lxi = toneLevel[i] - 10.0 * log10(den);
        if (!isfinite(lxi)) return -1;
        /* Masked tones contribute zero; nonfinite calculations are errors. */
        outLx[i] = lxi > 0.0 ? lxi : 0.0;
    }
    return 0;
}

static double tonal_weighting(const double *toneFreq, const double *toneBw,
                              const double *toneLx, const double *footprint, int nTones)
{
    double sumSq = 0.0;

    for (int i = 0; i < nTones; i++) {
        double fc = toneFreq[i];
        double bw = toneBw[i];
        double deltaL = toneLx[i];
        if (!isfinite(fc) || !isfinite(bw) || !isfinite(deltaL)
            || !isfinite(footprint[i]) || footprint[i] < 0.0) return NAN;
        if (fc <= 0.0 || bw <= 0.0 || deltaL <= 0.0) continue;

        double intrinsic = sqrt(fmax(bw*bw - footprint[i]*footprint[i], 0.0));
        double zup = fq_to_bark(fc + 0.5 * intrinsic);
        double zlow = fq_to_bark(fc - 0.5 * intrinsic);
        double dz = zup - zlow;
        double w1 = 0.13 / (dz + 0.13);
        double x = fc / 700.0 + 700.0 / fc;
        double w2 = pow(1.0 / sqrt(1.0 + 0.2 * x * x), 0.29);
        double w3 = pow(1.0 - exp(-deltaL / 15.0), 0.29);

        if (!isfinite(w1) || !isfinite(w2) || !isfinite(w3)
            || w1 < 0.0 || w2 < 0.0 || w3 < 0.0) return NAN;

        double ww1 = pow(w1, 1.0 / 0.29);
        double ww2 = pow(w2, 1.0 / 0.29);
        double ww3 = pow(w3, 1.0 / 0.29);
        double prod = ww1 * ww2 * ww3;
        sumSq += prod * prod;
    }

    return sqrt(sumSq);
}

#include "tonality_narrowband.h"

int tonality_aures1985(
    const double *signal,
    int           numSamples,
    double        sampleRate,
    int           soundField,
    double        timeSkip,
    double       *outTonality,
    double       *outTonalWeighting,
    double       *outLoudnessWeighting,
    double       *outTime,
    int          *pNumFrames,
    double       *outStats)
{
    const double *audio = signal;
    int audioLen = numSamples;
    double fs = sampleRate;

    if (!signal || !outTonality || !pNumFrames || numSamples <= 0) {
        return TonalityAures1985ErrorInvalidArgument;
    }
    if (sampleRate <= 0.0 || !isfinite(sampleRate)) {
        return TonalityAures1985ErrorInvalidSampleRate;
    }
    if (soundField != MS_SOUND_FIELD_FREE && soundField != MS_SOUND_FIELD_DIFFUSE) {
        return TonalityAures1985ErrorInvalidArgument;
    }

    if (fs != 48000.0) return TonalityAures1985ErrorInvalidSampleRate;

    int n = round_to_int(fs * AURES_TIME_RESOLUTION);
    int hop = round_to_int(0.5 * (double)n);
    int numFrames = audioLen < n ? 0 : 1 + (audioLen - n) / hop;
    if (n <= 0 || hop <= 0 || numFrames <= 0) {
        return TonalityAures1985ErrorSignalTooShort;
    }
    if (*pNumFrames < numFrames) {
        *pNumFrames = numFrames;
        return TonalityAures1985ErrorOutputTooSmall;
    }
    *pNumFrames = numFrames;

    cfft_plan fftPlan = make_cfft_plan((size_t)n);
    if (!fftPlan) {
        return TonalityAures1985ErrorFFTPlan;
    }

    double *window = (double*)calloc((size_t)n, sizeof(double));
    double *fftBuf = (double*)calloc((size_t)2 * n, sizeof(double));
    double *filtBuf = (double*)calloc((size_t)2 * n, sizeof(double));
    double *spl = (double*)calloc((size_t)n, sizeof(double));
    double *frame = (double*)calloc((size_t)n, sizeof(double));
    double *filtered = (double*)calloc((size_t)n, sizeof(double));

    int minIdx = (int)ceil(1.0 + AURES_MIN_FREQUENCY * ((double)n / fs)) - 1;
    int maxIdx = (int)ceil(1.0 + AURES_MAX_FREQUENCY * ((double)n / fs)) - 1;
    if (minIdx < 0) minIdx = 0;
    if (maxIdx >= n / 2) maxIdx = n / 2;
    int cropLen = maxIdx - minIdx + 1;

    int *toneBin = (int*)calloc((size_t)cropLen, sizeof(int));
    int *toneCropIdx = (int*)calloc((size_t)cropLen, sizeof(int));
    double *toneFreq = (double*)calloc((size_t)cropLen, sizeof(double));
    double *toneLevel = (double*)calloc((size_t)cropLen, sizeof(double));
    double *toneBw = (double*)calloc((size_t)cropLen, sizeof(double));
    double *toneLx = (double*)calloc((size_t)cropLen, sizeof(double));
    double *footprint = (double*)calloc((size_t)cropLen, sizeof(double));
    double *notchWidth = (double*)calloc((size_t)cropLen, sizeof(double));
    double *nbScratch = (double*)calloc((size_t)7*(size_t)(n/2+1), sizeof(double));
    unsigned char *taken = (unsigned char*)calloc((size_t)(n/2+1), 1u);

    if (!window || !fftBuf || !filtBuf || !spl || !frame || !filtered ||
        !toneBin || !toneCropIdx || !toneFreq || !toneLevel || !toneBw || !toneLx ||
        !footprint || !notchWidth || !nbScratch || !taken || cropLen < 8) {
        free(window); free(fftBuf); free(filtBuf); free(spl); free(frame); free(filtered);
        free(toneBin); free(toneCropIdx); free(toneFreq); free(toneLevel); free(toneBw); free(toneLx);
        free(footprint); free(notchWidth); free(nbScratch); free(taken);
        destroy_cfft_plan(fftPlan);
        return TonalityAures1985ErrorMemoryAlloc;
    }

    int result = numFrames;
    hann_window(window, n);
    double winMean = mean_array(window, n);
    double fftGain = sqrt(2.0) / ((double)n * winMean);
    double df = fs / (double)n;
    double windowPower = 0.0;
    for (int i = 0; i < n; ++i) windowPower += window[i]*window[i];
    double enbw = windowPower / ((double)n * winMean * winMean);
    int nearGuard = round_to_int(25.0/df), farGuard = round_to_int(37.5/df);
    if (nearGuard < 1) nearGuard = 1;
    if (farGuard < 2) farGuard = 2;
    double loudnessTimeSkip = AURES_TIME_RESOLUTION * AURES_LOUDNESS_SKIP_FRACTION;

    for (int fr = 0; fr < numFrames; fr++) {
        int start = fr * hop;
        int nTones = 0;
        aures_region regions[AURES_MAX_REGIONS];
        int nRegions;

        for (int i = 0; i < n; i++) {
            frame[i] = audio[start + i] * window[i];
            fftBuf[2 * i] = frame[i] * fftGain;
            fftBuf[2 * i + 1] = 0.0;
        }

        if (do_fft(fftBuf, n, fftPlan) != 0) { result = TonalityAures1985ErrorMemoryAlloc; goto cleanup; }
        for (int i = 0; i < n; i++) {
            double re = fftBuf[2 * i];
            double im = fftBuf[2 * i + 1];
            double energy = re * re + im * im;
            double relativeEnergy = (energy + AURES_TINY_VALUE) / AURES_I_REF;
            if (!isfinite(relativeEnergy)) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
            spl[i] = 10.0 * log10(relativeEnergy);
        }

        nRegions = find_narrowband(spl, n/2+1, df, nbScratch, taken, regions);
        if (nRegions < 0) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
        for (int ci = 0; ci < cropLen; ci++) {
            int bi = minIdx + ci;
            if (bi-farGuard < 0 || bi+farGuard >= n) continue;
            double p = spl[bi];
            if (p > spl[bi - 1] &&
                p >= spl[bi + 1] &&
                p - spl[bi - farGuard] >= AURES_TONE_THRESHOLD &&
                p - spl[bi - nearGuard] >= AURES_TONE_THRESHOLD &&
                p - spl[bi + nearGuard] >= AURES_TONE_THRESHOLD &&
                p - spl[bi + farGuard] >= AURES_TONE_THRESHOLD) {
                int absorbed = 0;
                for (int r = 0; r < nRegions; ++r)
                    if (bi >= regions[r].first && bi <= regions[r].last) absorbed = 1;
                if (absorbed || p <= 0.0) continue;
                toneBin[nTones] = bi;
                toneCropIdx[nTones] = ci;
                toneLevel[nTones] = p;
                toneFreq[nTones] = (double)bi * fs / (double)n;
                nTones++;
            }
        }

        if (nTones > 0) {
            if (estimate_bandwidths(spl, n, fs, toneBin, toneLevel, nTones, toneBw) != 0) {
                result = TonalityAures1985ErrorNumerical;
                goto cleanup;
            }
        }

        double wTonal = 0.0;
        double wLoudness = 0.0;
        double tonality = 0.0;

        if (nTones > 0 || nRegions > 0) {
            for (int i = 0; i < n; i++) {
                filtBuf[2 * i] = frame[i];
                filtBuf[2 * i + 1] = 0.0;
            }
            if (do_fft(filtBuf, n, fftPlan) != 0) { result = TonalityAures1985ErrorMemoryAlloc; goto cleanup; }
            for (int i = 0; i < nTones; ++i) {
                footprint[i] = 1.43*df; /* SQAT worst-case sinusoid estimator width. */
                notchWidth[i] = fmax(toneBw[i], 4.0*df); /* Hann main-lobe coverage. */
            }
            remove_tones_preserve_phase(filtBuf, n, fs, toneFreq, notchWidth, nTones);
            for (int r = 0; r < nRegions; ++r) {
                for (int k = regions[r].first; k <= regions[r].last; ++k) {
                    double magnitude = sqrt(pow(10.0, region_bridge(&regions[r], k)/10.0)
                                           * AURES_I_REF)/fftGain;
                    double phase = atan2(filtBuf[2*k+1], filtBuf[2*k]);
                    filtBuf[2*k] = magnitude*cos(phase);
                    filtBuf[2*k+1] = magnitude*sin(phase);
                }
                if (nTones >= cropLen) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
                toneFreq[nTones] = regions[r].frequency;
                toneBw[nTones] = regions[r].width;
                toneLevel[nTones] = regions[r].level - 10.0*log10(enbw);
                footprint[nTones++] = 2.667*df; /* SQAT region estimator footprint. */
            }
            /* Stable frequency ordering determines the secondary masking slopes. */
            for (int i = 1; i < nTones; ++i) {
                for (int j = i; j > 0 && toneFreq[j] < toneFreq[j-1]; --j) {
                    double swap;
                    swap = toneFreq[j]; toneFreq[j] = toneFreq[j-1]; toneFreq[j-1] = swap;
                    swap = toneBw[j]; toneBw[j] = toneBw[j-1]; toneBw[j-1] = swap;
                    swap = toneLevel[j]; toneLevel[j] = toneLevel[j-1]; toneLevel[j-1] = swap;
                    swap = footprint[j]; footprint[j] = footprint[j-1]; footprint[j-1] = swap;
                }
            }
            for (int i = 0; i < nTones; ++i) {
                toneCropIdx[i] = -1; /* No five-bin exclusion unless centre is on this grid. */
                double z = fq_to_bark(toneFreq[i]);
                for (int k = 0; k < cropLen; ++k)
                    if (fq_to_bark((minIdx+k)*df) == z) { toneCropIdx[i] = k; break; }
            }
            /* Residual spectrum supplies noise masking; then restore symmetry. */
            for (int k = 0; k <= n/2; ++k) {
                double magnitude = complex_abs_at(filtBuf, k)*fftGain;
                spl[k] = 10.0*log10(magnitude*magnitude/AURES_I_REF + AURES_TINY_VALUE);
                if (!isfinite(spl[k])) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
            }
            for (int k = 1; k < n/2; ++k) {
                filtBuf[2*(n-k)] = filtBuf[2*k];
                filtBuf[2*(n-k)+1] = -filtBuf[2*k+1];
            }
            if (do_ifft(filtBuf, n, fftPlan) != 0) { result = TonalityAures1985ErrorMemoryAlloc; goto cleanup; }
            for (int i = 0; i < n; i++) filtered[i] = filtBuf[2 * i];

            double lTotal = 0.0;
            double lFiltered = 0.0;
            int retTotal = compute_stationary_loudness(frame, n, fs, soundField,
                                                       loudnessTimeSkip, &lTotal);
            int retFiltered = compute_stationary_loudness(filtered, n, fs, soundField,
                                                          loudnessTimeSkip, &lFiltered);
            if (retTotal < 0 || retFiltered < 0) {
                result = TonalityAures1985ErrorLoudnessFailed;
                goto cleanup;
            }

            if (lTotal > AURES_TINY_VALUE) {
                wLoudness = 1.0 - (lFiltered / lTotal);
                if (!isfinite(wLoudness)) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
                /* Spectral replacement can increase residual loudness. The
                 * reference then defines the tonal contribution as zero. */
                if (wLoudness < 0.0) wLoudness = 0.0;
            }

            if (spl_excess(spl + minIdx, cropLen, minIdx, fs, n,
                          toneFreq, toneLevel, toneCropIdx, nTones, toneLx) != 0) {
                result = TonalityAures1985ErrorNumerical;
                goto cleanup;
            }
            wTonal = tonal_weighting(toneFreq, toneBw, toneLx, footprint, nTones);
            tonality = AURES_REFERENCE_SCALE * pow(wTonal, 0.29) * pow(wLoudness, 0.79);
            if (!isfinite(tonality)) { result = TonalityAures1985ErrorNumerical; goto cleanup; }
        }

        outTonality[fr] = tonality;
        if (outTonalWeighting) outTonalWeighting[fr] = wTonal;
        if (outLoudnessWeighting) outLoudnessWeighting[fr] = wLoudness;
        if (outTime) outTime[fr] = ((double)start + 1.0) / fs;
    }

    int statsRet = 0;
    if (outStats) {
        if (outTime) {
            statsRet = fill_stats(outTonality, outTime, numFrames, timeSkip, outStats);
        } else {
            double *tmpTime = (double*)calloc((size_t)numFrames, sizeof(double));
            if (!tmpTime) {
                statsRet = TonalityAures1985ErrorMemoryAlloc;
            } else {
                for (int fr = 0; fr < numFrames; fr++) {
                    tmpTime[fr] = ((double)(fr * hop) + 1.0) / fs;
                }
                statsRet = fill_stats(outTonality, tmpTime, numFrames, timeSkip, outStats);
                free(tmpTime);
            }
        }
    }

    if (statsRet < 0) result = statsRet;
cleanup:
    free(window);
    free(fftBuf);
    free(filtBuf);
    free(spl);
    free(frame);
    free(filtered);
    free(toneBin);
    free(toneCropIdx);
    free(toneFreq);
    free(toneLevel);
    free(toneBw);
    free(toneLx);
    free(footprint); free(notchWidth); free(nbScratch); free(taken);
    destroy_cfft_plan(fftPlan);

    return result;
}
