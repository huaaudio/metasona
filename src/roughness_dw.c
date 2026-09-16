// SPDX-License-Identifier: Apache-2.0 AND MIT
// MetaSona adaptation: Jiahua Zhang and Codex, September 2026.
// Based on the user-supplied C port; aligned to SQAT e6228b789fc9.
// Model sources: Roughness_Daniel1997, Get_Hweight_roughness,
// Get_gzi_roughness, Terhardt_filterbank and Terhardt_filterbank_params.
// MetaSona changes: complete frames, initialized cleanup and checked FFT/numeric errors.
// Model/data provenance: NOTICE, THIRD_PARTY.md, LICENSES/MIT-SQAT.txt.

/************************************************************************/
/*  Roughness calculation according to Daniel & Weber (1997)            */
/*  "Psychoacoustical roughness: implementation of an optimized model"  */
/*  Acustica(83), 113-123.                                              */
/*                                                                      */
/*  Ported from MATLAB (Roughness_Daniel1997.m, SQAT toolbox)           */
/************************************************************************/

#include "roughness_dw.h"
#include "pocketfft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ======================================================================
   Internal constants
   ====================================================================== */

#define TIME_RESOLUTION  0.2   /* seconds per frame */
#define N_CHNO           47    /* number of critical band channels */
/* SQAT's documented 48 kHz reference calibration, not a fit to our tests. */
#define SQAT_REFERENCE_ASPER 1.006197
#define DZ 0.5
#define CAL (0.25 / SQAT_REFERENCE_ASPER / DZ)

/* Bark table: [BarkNum, Fc_low, Fc_center, BarkNum+0.5] - 25 rows */
static const double BarkTable[25][4] = {
    { 0,     0,    50,  0.5},
    { 1,   100,   150,  1.5},
    { 2,   200,   250,  2.5},
    { 3,   300,   350,  3.5},
    { 4,   400,   450,  4.5},
    { 5,   510,   570,  5.5},
    { 6,   630,   700,  6.5},
    { 7,   770,   840,  7.5},
    { 8,   920,  1000,  8.5},
    { 9,  1080,  1170,  9.5},
    {10,  1270,  1370, 10.5},
    {11,  1480,  1600, 11.5},
    {12,  1720,  1850, 12.5},
    {13,  2000,  2150, 13.5},
    {14,  2320,  2500, 14.5},
    {15,  2700,  2900, 15.5},
    {16,  3150,  3400, 16.5},
    {17,  3700,  4000, 17.5},
    {18,  4400,  4800, 18.5},
    {19,  5300,  5800, 19.5},
    {20,  6400,  7000, 20.5},
    {21,  7700,  8500, 21.5},
    {22,  9500, 10500, 22.5},
    {23, 12000, 13500, 23.5},
    {24, 15500, 20000, 24.5}
};

/* Hearing threshold table */
static const double HTres[][2] = {
    { 0,      130}, { 0.01,    70}, { 0.17,    60}, { 0.8,     30},
    { 1,       25}, { 1.5,     20}, { 2,       15}, { 3.3,     10},
    { 4,      8.1}, { 5,      6.3}, { 6,        5}, { 8,      3.5},
    {10,      2.5}, {12,      1.7}, {13.3,      0}, {15,     -2.5},
    {16,       -4}, {17,     -3.7}, {18,     -1.5}, {19,      1.4},
    {20,      3.8}, {21,        5}, {22,      7.5}, {23,       15},
    {24,       48}, {24.5,     60}, {25,      130}
};
#define N_HTRES  27

/* a0 table */
static const double a0tab[][2] = {
    { 0,       0}, {10,       0}, {12,    1.15}, {13,    2.31},
    {14,    3.85}, {15,    5.62}, {16,    6.92}, {16.5,  7.38},
    {17,    6.92}, {18,    4.23}, {18.5,  2.31}, {19,       0},
    {20,   -1.43}, {21,   -2.59}, {21.5, -3.57}, {22,   -5.19},
    {22.5, -7.41}, {23,   -11.3}, {23.5,   -20}, {24,     -40},
    {25,    -130}, {26,    -999}
};
#define N_A0TAB  22

/* gr table for gzi calculation (Aures modulation depth weighting, from MoSQITo) */
static const double gr_x[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24};
static const double gr_y[] = {0.15,0.26,0.38,0.47,0.54,0.65,0.76,0.83,0.90,0.98,0.98,0.90,0.80,0.70,0.62,0.54,0.49,0.43,0.39,0.35,0.30,0.30,0.30,0.30,0.30};
#define N_GR  25

/* H-weight tables */
static const double H2_tab[][2] = {
    {0,0},{17,0.8},{23,0.95},{25,0.975},{32,1},{37,0.975},{48,0.9},
    {67,0.8},{90,0.7},{114,0.6},{171,0.4},{206,0.3},{247,0.2},{294,0.1},{358,0}
};
#define N_H2 15

static const double H5_tab[][2] = {
    {0,0},{32,0.8},{43,0.95},{56,1},{69,0.975},{92,0.9},{120,0.8},
    {142,0.7},{165,0.6},{231,0.4},{277,0.3},{331,0.2},{397,0.1},{502,0}
};
#define N_H5 14

static const double H16_tab[][2] = {
    {0,0},{23.5,0.4},{34,0.6},{47,0.8},{56,0.9},{63,0.95},{79,1},
    {100,0.975},{115,0.95},{135,0.9},{159,0.85},{172,0.8},{194,0.7},
    {215,0.6},{244,0.5},{290,0.4},{348,0.3},{415,0.2},{500,0.1},{645,0}
};
#define N_H16 20

static const double H21_tab[][2] = {
    {0,0},{19,0.4},{44,0.8},{52.5,0.9},{58,0.95},{75,1},{101.5,0.95},
    {114.5,0.9},{132.5,0.85},{143.5,0.8},{165.5,0.7},{197.5,0.6},
    {241,0.5},{290,0.4},{348,0.3},{415,0.2},{500,0.1},{645,0}
};
#define N_H21 18

static const double H42_tab[][2] = {
    {0,0},{15,0.4},{41,0.8},{49,0.9},{53,0.965},{64,0.99},{71,1},
    {88,0.95},{94,0.9},{106,0.85},{115,0.8},{137,0.7},{180,0.6},
    {238,0.5},{290,0.4},{348,0.3},{415,0.2},{500,0.1},{645,0}
};
#define N_H42 19

/* ======================================================================
   Utility functions
   ====================================================================== */

static double interp1_table(const double table[][2], int n, double x)
{
    if (x <= table[0][0]) return table[0][1];
    if (x >= table[n-1][0]) return table[n-1][1];
    for (int i = 0; i < n - 1; i++) {
        if (x >= table[i][0] && x <= table[i+1][0]) {
            double t = (x - table[i][0]) / (table[i+1][0] - table[i][0]);
            return table[i][1] + t * (table[i+1][1] - table[i][1]);
        }
    }
    return table[n-1][1];
}

static double interp1_xy(const double *xarr, const double *yarr, int n, double x)
{
    if (x <= xarr[0]) return yarr[0];
    if (x >= xarr[n-1]) return yarr[n-1];
    for (int i = 0; i < n - 1; i++) {
        if (x >= xarr[i] && x <= xarr[i+1]) {
            double t = (x - xarr[i]) / (xarr[i+1] - xarr[i]);
            return yarr[i] + t * (yarr[i+1] - yarr[i]);
        }
    }
    return yarr[n-1];
}

static double db2mag(double dB) { return pow(10.0, dB / 20.0); }
static double mag2db(double mag) { return (mag <= 0.0) ? -400.0 : 20.0 * log10(mag); }

static void blackman_window(double *w, int N)
{
    for (int i = 0; i < N; i++)
        w[i] = 0.42 - 0.5*cos(2.0*M_PI*i/N) + 0.08*cos(4.0*M_PI*i/N);
}

static double arr_mean(const double *a, int n)
{
    double s = 0; for (int i = 0; i < n; i++) s += a[i]; return s / n;
}

/* ======================================================================
   FFT wrappers using pocketfft (cfft_plan).
   Data: interleaved [re0,im0,re1,im1,...] of length 2*N.
   ====================================================================== */

/* Forward FFT: no scaling (fct=1.0) */
static int do_fft(double *data, int N, cfft_plan plan)
{
    (void)N; /* Transform length is owned by the PocketFFT plan. */
    return cfft_forward(plan, data, 1.0);
}

/* Inverse FFT: result scaled by 1/N */
static int do_ifft(double *data, int N, cfft_plan plan)
{
    return cfft_backward(plan, data, 1.0/(double)N);
}

/* ======================================================================
   Simple linear resampling
   ====================================================================== */
static double* resample_linear(const double *in, int inLen, int inRate,
                               int outRate, int *outLen)
{
    double ratio = (double)inRate / (double)outRate;
    *outLen = (int)((double)inLen / ratio);
    double *out = (double*)calloc(*outLen, sizeof(double));
    if (!out) return NULL;
    for (int i = 0; i < *outLen; i++) {
        double s = i * ratio;
        int idx = (int)s;
        double f = s - idx;
        if (idx+1 < inLen) out[i] = in[idx]*(1.0-f) + in[idx+1]*f;
        else if (idx < inLen) out[i] = in[idx];
    }
    return out;
}

/* ======================================================================
   Main roughness computation
   ====================================================================== */

int roughness_dw(
    const double *pSignal,
    int           numSamples,
    double        sampleRate,
    double       *outRoughness,
    double       *outSpecRoughness,
    int          *pNumFrames)
{
    double *Bark2_freq = NULL;
    double *Bark2_bark = NULL;
    double *Barkno = NULL;
    int *zb = NULL;
    double *MinExcdB = NULL;
    double *a0 = NULL;
    double *Hweight = NULL;
    double *fft_buf = NULL;
    double *hBPi_all = NULL;
    double *dataIn = NULL;
    double *TempIn_re = NULL;
    double *TempIn_im = NULL;
    double *Slopes = NULL;
    int *whichL = NULL;
    double *Lg = NULL;
    double *LdB = NULL;
    double *S2 = NULL;
    int *wzF = NULL;
    int *wzC = NULL;
    double **t_ifft_buf = NULL;
    double **t_ei_k = NULL;
    double **t_etmp_abs = NULL;
    double **t_fft_etmp = NULL;
    cfft_plan *t_fft_plan = NULL;
    int nThreads = 1;
    int computationError = RoughnessErrorMemoryAlloc;
    const double *audio = pSignal;
    int audioLen = numSamples;
    double fs = sampleRate;
    double *resampledAudio = NULL;

    /* Resample to 48kHz if not a supported rate */
    if (fs != 44100.0 && fs != 40960.0 && fs != 48000.0) {
        int newLen;
        resampledAudio = resample_linear(pSignal, numSamples, (int)fs, 48000, &newLen);
        if (!resampledAudio) return RoughnessErrorMemoryAlloc;
        audio = resampledAudio;
        audioLen = newLen;
        fs = 48000.0;
    }

    /* N = window length = FFT length (matches MATLAB exactly) */
    int N = (int)(fs * TIME_RESOLUTION);
    int hopsize = N / 2;
    int numFrames = audioLen < N ? 0 : 1 + (audioLen - N) / hopsize;

    if (numFrames <= 0) {
        if (resampledAudio) free(resampledAudio);
        return RoughnessErrorSignalTooShort;
    }
    if (*pNumFrames < numFrames) {
        if (resampledAudio) free(resampledAudio);
        return RoughnessErrorOutputTooSmall;
    }
    *pNumFrames = numFrames;

    /* Create pocketfft plan */
    cfft_plan fft_plan = make_cfft_plan((size_t)N);
    if (!fft_plan) {
        if (resampledAudio) free(resampledAudio);
        return RoughnessErrorMemoryAlloc;
    }

    /* Precompute window */
    double *window = (double*)calloc(N, sizeof(double));
    if (!window) { destroy_cfft_plan(fft_plan); if (resampledAudio) free(resampledAudio); return RoughnessErrorMemoryAlloc; }
    blackman_window(window, N);
    double winMean = arr_mean(window, N);
    double AmpCal = (1.0 / (20e-6 * sqrt(2.0))) * 2.0 / (N * winMean);

    /* Frequency parameters */
    int N2 = N / 2 + 1;
    double dFs = fs / N;
    int N0 = (int)round(20.0 * N / fs);
    int Ntop = (int)round(20000.0 * N / fs);
    if (Ntop >= N2) Ntop = N2 - 1;

    /* Build Bark2 sorted table */
    int nBark2 = 50;
    Bark2_freq = (double*)calloc(nBark2, sizeof(double));
    Bark2_bark = (double*)calloc(nBark2, sizeof(double));
    if (!Bark2_freq || !Bark2_bark) goto fail_alloc;

    for (int i = 0; i < 25; i++) {
        Bark2_freq[2*i]   = BarkTable[i][1];
        Bark2_bark[2*i]   = BarkTable[i][0];
        Bark2_freq[2*i+1] = BarkTable[i][2];
        Bark2_bark[2*i+1] = BarkTable[i][3];
    }
    /* Sort by frequency (bubble sort, only 50 elements) */
    for (int i = 0; i < nBark2-1; i++)
        for (int j = i+1; j < nBark2; j++)
            if (Bark2_freq[j] < Bark2_freq[i]) {
                double t = Bark2_freq[i]; Bark2_freq[i] = Bark2_freq[j]; Bark2_freq[j] = t;
                t = Bark2_bark[i]; Bark2_bark[i] = Bark2_bark[j]; Bark2_bark[j] = t;
            }

    /* Barkno array */
    Barkno = (double*)calloc(N2, sizeof(double));
    if (!Barkno) goto fail_alloc;
    for (int f = N0; f <= Ntop; f++)
        Barkno[f] = interp1_xy(Bark2_freq, Bark2_bark, nBark2, f*dFs);

    /* Cf and Bf */
    int Cf[24], Bf_start[25];
    for (int a = 0; a < 24; a++)
        Cf[a] = (int)round(BarkTable[a+1][1] * N / fs) - N0;
    Bf_start[0] = (int)round(BarkTable[0][2] * N / fs) - N0;
    for (int a = 0; a < 24; a++)
        Bf_start[a+1] = (int)round(BarkTable[a+1][2] * N / fs) - N0;

    /* zb sorted */
    int nzb = 49;
    zb = (int*)calloc(nzb, sizeof(int));
    if (!zb) goto fail_alloc;
    for (int i = 0; i < 25; i++) zb[i] = Bf_start[i];
    for (int i = 0; i < 24; i++) zb[25+i] = Cf[i];
    for (int i = 0; i < nzb-1; i++)
        for (int j = i+1; j < nzb; j++)
            if (zb[j] < zb[i]) { int t = zb[i]; zb[i] = zb[j]; zb[j] = t; }

    /* MinExcdB */
    int nBins = Ntop - N0 + 1;
    MinExcdB = (double*)calloc(nBins, sizeof(double));
    if (!MinExcdB) goto fail_alloc;
    for (int i = 0; i < nBins; i++)
        MinExcdB[i] = interp1_table(HTres, N_HTRES, Barkno[N0+i]);

    /* SQAT samples the threshold at its tabulated channel frequencies. */
    double MinBf[N_CHNO];
    for (int i = 0; i < N_CHNO; i++)
        MinBf[i] = MinExcdB[zb[i]];

    /* a0 */
    a0 = (double*)calloc(N, sizeof(double));
    if (!a0) goto fail_alloc;
    for (int i = 0; i < N; i++) a0[i] = 1.0;
    for (int k = N0; k <= Ntop; k++)
        a0[k] = db2mag(interp1_table(a0tab, N_A0TAB, Barkno[k]));

    /* SQAT squares sqrt(g)*m*k; multiplying g after squaring is equivalent. */
    double gzi[N_CHNO];
    for (int k = 0; k < N_CHNO; k++)
        gzi[k] = interp1_xy(gr_x, gr_y, N_GR, (k+1)/2.0);

    /* H-weight matrix [N_CHNO x N] */
    int DCbins = 2;
    int lastH2  = (int)floor(358.0/fs*N);
    int lastH5  = (int)floor(502.0/fs*N);
    int lastH16 = (int)floor(645.0/fs*N);

    Hweight = (double*)calloc((size_t)N_CHNO * N, sizeof(double));
    if (!Hweight) goto fail_alloc;
    #define HW(ch,bin) Hweight[(ch)*N + (bin)]

    for (int k = DCbins; k < lastH2 && k < N; k++)
        HW(1,k) = interp1_table(H2_tab, N_H2, (double)k*fs/N);
    for (int k = DCbins; k < lastH5 && k < N; k++)
        HW(4,k) = interp1_table(H5_tab, N_H5, (double)k*fs/N);
    for (int k = DCbins; k < lastH16 && k < N; k++)
        HW(15,k) = interp1_table(H16_tab, N_H16, (double)k*fs/N);
    for (int k = DCbins; k < lastH16 && k < N; k++)
        HW(20,k) = interp1_table(H21_tab, N_H21, (double)k*fs/N);
    for (int k = DCbins; k < lastH16 && k < N; k++)
        HW(41,k) = interp1_table(H42_tab, N_H42, (double)k*fs/N);

    /* Copy channel weights */
    for (int k = 0; k < N; k++) { HW(0,k)=HW(1,k); HW(2,k)=HW(1,k); HW(3,k)=HW(1,k); }
    for (int l = 5; l <= 14; l++) for (int k = 0; k < N; k++) HW(l,k) = HW(4,k);
    for (int l = 16; l <= 19; l++) for (int k = 0; k < N; k++) HW(l,k) = HW(15,k);
    for (int l = 21; l <= 40; l++) for (int k = 0; k < N; k++) HW(l,k) = HW(20,k);
    for (int l = 42; l <= 46; l++) for (int k = 0; k < N; k++) HW(l,k) = HW(41,k);

    /* Working arrays (per-frame, shared across channels) */
    fft_buf = (double*)calloc(2*N, sizeof(double));
    hBPi_all = (double*)calloc((size_t)N_CHNO * N, sizeof(double));
    dataIn = (double*)calloc(N, sizeof(double));
    TempIn_re = (double*)calloc(N, sizeof(double));
    TempIn_im = (double*)calloc(N, sizeof(double));
    Slopes = (double*)calloc((size_t)nBins * N_CHNO, sizeof(double));

    /* Per-frame arrays hoisted out of the loop (max size = nBins) */
    whichL = (int*)calloc(nBins, sizeof(int));
    Lg = (double*)calloc(nBins, sizeof(double));
    LdB = (double*)calloc(nBins, sizeof(double));
    S2 = (double*)calloc(nBins > 0 ? nBins : 1, sizeof(double));
    wzF = (int*)calloc(nBins > 0 ? nBins : 1, sizeof(int));
    wzC = (int*)calloc(nBins > 0 ? nBins : 1, sizeof(int));

    /* Per-thread work buffers for the parallel channel loop */
#ifdef _OPENMP
    nThreads = omp_get_max_threads();
#endif
    t_ifft_buf = (double**)calloc(nThreads, sizeof(double*));
    t_ei_k = (double**)calloc(nThreads, sizeof(double*));
    t_etmp_abs = (double**)calloc(nThreads, sizeof(double*));
    t_fft_etmp = (double**)calloc(nThreads, sizeof(double*));
    t_fft_plan = (cfft_plan*)calloc(nThreads, sizeof(cfft_plan));
    int thread_alloc_ok = (t_ifft_buf && t_ei_k && t_etmp_abs && t_fft_etmp && t_fft_plan);
    if (thread_alloc_ok) {
        for (int t = 0; t < nThreads; t++) {
            t_ifft_buf[t] = (double*)calloc(2*N, sizeof(double));
            t_ei_k[t]     = (double*)calloc(N, sizeof(double));
            t_etmp_abs[t] = (double*)calloc(N, sizeof(double));
            t_fft_etmp[t] = (double*)calloc(2*N, sizeof(double));
            t_fft_plan[t] = make_cfft_plan((size_t)N);
            if (!t_ifft_buf[t] || !t_ei_k[t] || !t_etmp_abs[t] || !t_fft_etmp[t] ||
                !t_fft_plan[t]) {
                thread_alloc_ok = 0; break;
            }
        }
    }

    if (!fft_buf||!hBPi_all||
        !dataIn||!TempIn_re||!TempIn_im||!Slopes||
        !whichL||!Lg||!LdB||!S2||!wzF||!wzC||!thread_alloc_ok)
        goto fail_alloc;

    #define HBPI(ch,bin) hBPi_all[(ch)*N + (bin)]
    #define SLP(w,ch)    Slopes[(w)*N_CHNO + (ch)]

    double specR_sum[N_CHNO];
    memset(specR_sum, 0, sizeof(specR_sum));

    /* ===== Main processing loop ===== */
    for (int fr = 0; fr < numFrames; fr++) {
        int startIdx = fr * hopsize;

        /* Window input */
        for (int i = 0; i < N; i++)
            dataIn[i] = audio[startIdx + i] * window[i];

        /* FFT with AmpCal scaling and a0 weighting */
        for (int i = 0; i < N; i++) { fft_buf[2*i] = dataIn[i]*AmpCal; fft_buf[2*i+1] = 0; }
        if (do_fft(fft_buf, N, fft_plan) != 0) goto fail_alloc;
        for (int i = 0; i < N; i++) {
            TempIn_re[i] = fft_buf[2*i]   * a0[i];
            TempIn_im[i] = fft_buf[2*i+1] * a0[i];
        }

        /* Magnitude spectrum in hearing range */
        int sizL = 0;

        for (int i = 0; i < nBins; i++) {
            int k = N0 + i;
            double re = TempIn_re[k], im = TempIn_im[k];
            Lg[i] = sqrt(re*re + im*im);
            if (!isfinite(Lg[i])) {
                computationError = RoughnessErrorNumerical;
                goto fail_alloc;
            }
            LdB[i] = mag2db(Lg[i]);
            if (LdB[i] > MinExcdB[i]) whichL[sizL++] = i;
        }

        /* Slopes */
        memset(Slopes, 0, (size_t)nBins * N_CHNO * sizeof(double));
        double S1 = -27.0;

        for (int w = 0; w < sizL; w++) {
            int idx = whichL[w];
            double freq_w = (double)(N0+idx) * fs / N;
            double steep = -24.0 - 230.0/freq_w + 0.2*LdB[idx];
            S2[w] = (steep < 0) ? steep : 0.0;
            int bi = N0 + idx;
            if (bi >= 0 && bi < N2) {
                wzF[w] = (int)floor(2.0 * Barkno[bi]);
                wzC[w] = (int)ceil(2.0 * Barkno[bi]);
            }
        }

        for (int w = 0; w < sizL; w++) {
            int idx = whichL[w];
            double Ltmp = LdB[idx];
            int bi = N0 + idx;
            double Btmp = (bi >= 0 && bi < N2) ? Barkno[bi] : 0.0;

            for (int l = 0; l < wzF[w] && l < N_CHNO; l++) {
                double St = S1*(Btmp - (l+1)*0.5) + Ltmp;
                if (St > MinBf[l]) SLP(w,l) = db2mag(St);
            }
            /* Upper slope starts at wzC-1 to match MoSQITo's ch_high = ceil(2*bark)-1 */
            {
                int upper_start = (wzC[w] > 0) ? wzC[w] - 1 : 0;
                for (int l = upper_start; l < N_CHNO; l++) {
                    double St = S2[w]*((l+1)*0.5 - Btmp) + Ltmp;
                    if (St > MinBf[l]) SLP(w,l) = db2mag(St);
                }
            }
        }

        int frameFftFailed = 0, frameNumericalFailed = 0;
        /* Per-channel processing */
        double h0[N_CHNO], hBPrms[N_CHNO], mdept[N_CHNO];
        memset(h0, 0, sizeof(h0));
        memset(hBPrms, 0, sizeof(hBPrms));
        memset(mdept, 0, sizeof(mdept));
        memset(hBPi_all, 0, (size_t)N_CHNO * N * sizeof(double));

        #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int ch = 0; ch < N_CHNO; ch++) {
            #ifdef _OPENMP
            int tid = omp_get_thread_num();
            #else
            int tid = 0;
            #endif
            double *my_ifft  = t_ifft_buf[tid];
            double *my_ei    = t_ei_k[tid];
            double *my_eabs  = t_etmp_abs[tid];
            double *my_fetmp = t_fft_etmp[tid];
            cfft_plan my_plan = t_fft_plan[tid];

            memset(my_ifft, 0, 2*N*sizeof(double));

            for (int w = 0; w < sizL; w++) {
                int idx = whichL[w];
                int fbin = N0 + idx;
                double ea = 0.0;

                /* SQAT's boundary channels have only their available
                 * neighbouring slope; all channel indices here are zero based. */
                if (ch > 0 && Lg[idx] > 0.0) ea = SLP(w, ch-1) / Lg[idx];
                if (ch < N_CHNO-1 && wzC[w] > ch+1 && Lg[idx] > 0.0)
                    ea = SLP(w, ch+1) / Lg[idx];
                if (wzF[w] == ch+1 || wzC[w] == ch+1) ea = 1.0;

                my_ifft[2*fbin]   += ea * TempIn_re[fbin];
                my_ifft[2*fbin+1] += ea * TempIn_im[fbin];
            }

            /* SQAT analytic half-spectrum: ei = 2*N*real(ifft(etmp)). */
            if (do_ifft(my_ifft, N, my_plan) != 0) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                frameFftFailed = 1;
                continue;
            }
            double sum_abs = 0;
            for (int i = 0; i < N; i++) {
                my_ei[i] = 2.0 * (double)N * my_ifft[2*i];
                my_eabs[i] = fabs(my_ei[i]);
                sum_abs += my_eabs[i];
            }
            h0[ch] = sum_abs / N;

            /* FFT of envelope, apply H-weight, IFFT */
            for (int i = 0; i < N; i++) { my_fetmp[2*i] = my_eabs[i]-h0[ch]; my_fetmp[2*i+1] = 0; }
            if (do_fft(my_fetmp, N, my_plan) != 0) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                frameFftFailed = 1;
                continue;
            }
            for (int i = 0; i < N; i++) { my_fetmp[2*i] *= HW(ch,i); my_fetmp[2*i+1] *= HW(ch,i); }
            if (do_ifft(my_fetmp, N, my_plan) != 0) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                frameFftFailed = 1;
                continue;
            }

            double rms_sq = 0;
            for (int i = 0; i < N; i++) {
                double v = 2.0 * my_fetmp[2*i];  /* do_ifft already divides by N */
                HBPI(ch,i) = v;
                rms_sq += v*v;
            }
            if (!isfinite(sum_abs) || !isfinite(rms_sq)) {
                #ifdef _OPENMP
                #pragma omp atomic write
                #endif
                frameNumericalFailed = 1;
                continue;
            }
            hBPrms[ch] = sqrt(rms_sq / N);

            if (h0[ch] > 0) {
                mdept[ch] = hBPrms[ch] / h0[ch];
                if (mdept[ch] > 1.0) mdept[ch] = 1.0;
            }
        }

        if (frameFftFailed) goto fail_alloc;
        if (frameNumericalFailed) {
            computationError = RoughnessErrorNumerical;
            goto fail_alloc;
        }
        /* Cross-correlation ki */
        double ki[45];
        for (int k = 0; k < 45; k++) {
            double m1 = 0, m2 = 0;
            for (int i = 0; i < N; i++) { m1 += HBPI(k,i); m2 += HBPI(k+2,i); }
            m1 /= N; m2 /= N;
            double c12 = 0, v1 = 0, v2 = 0;
            for (int i = 0; i < N; i++) {
                double d1 = HBPI(k,i)-m1, d2 = HBPI(k+2,i)-m2;
                c12 += d1*d2; v1 += d1*d1; v2 += d2*d2;
            }
            double den = sqrt(v1*v2);
            if (!isfinite(den)) {
                computationError = RoughnessErrorNumerical;
                goto fail_alloc;
            }
            ki[k] = (den > 0) ? c12/den : 0.0;
        }

        /* Specific roughness & total (MoSQITo: R_spec = gzi * (mdept * ki_product)^2) */
        double ri[N_CHNO];
        ri[0] = gzi[0] * (mdept[0]*ki[0]) * (mdept[0]*ki[0]);
        ri[1] = gzi[1] * (mdept[1]*ki[1]) * (mdept[1]*ki[1]);
        for (int k = 2; k < 45; k++) {
            double mki = mdept[k]*ki[k-2]*ki[k];
            ri[k] = gzi[k] * mki * mki;
        }
        ri[45] = gzi[45] * (mdept[45]*ki[43]) * (mdept[45]*ki[43]);
        ri[46] = gzi[46] * (mdept[46]*ki[44]) * (mdept[46]*ki[44]);

        double R = 0;
        for (int k = 0; k < N_CHNO; k++) {
            ri[k] *= CAL;
            R += ri[k];
            specR_sum[k] += ri[k];
        }
        if (!isfinite(R)) {
            computationError = RoughnessErrorNumerical;
            goto fail_alloc;
        }
        outRoughness[fr] = R * DZ;

    }

    /* Time-averaged specific roughness */
    if (outSpecRoughness) {
        for (int k = 0; k < N_CHNO; k++)
            outSpecRoughness[k] = specR_sum[k] / numFrames;
    }

    /* Cleanup */
    free(window); free(Bark2_freq); free(Bark2_bark); free(Barkno); free(zb);
    free(MinExcdB); free(a0); free(Hweight);
    free(fft_buf);
    free(hBPi_all); free(dataIn); free(TempIn_re); free(TempIn_im);
    free(Slopes);
    free(whichL); free(Lg); free(LdB);
    free(S2); free(wzF); free(wzC);
    for (int t = 0; t < nThreads; t++) {
        if (t_ifft_buf && t_ifft_buf[t]) free(t_ifft_buf[t]);
        if (t_ei_k && t_ei_k[t]) free(t_ei_k[t]);
        if (t_etmp_abs && t_etmp_abs[t]) free(t_etmp_abs[t]);
        if (t_fft_etmp && t_fft_etmp[t]) free(t_fft_etmp[t]);
        if (t_fft_plan && t_fft_plan[t]) destroy_cfft_plan(t_fft_plan[t]);
    }
    free(t_ifft_buf); free(t_ei_k); free(t_etmp_abs); free(t_fft_etmp); free(t_fft_plan);
    destroy_cfft_plan(fft_plan);
    if (resampledAudio) free(resampledAudio);
    #undef HW
    #undef HBPI
    #undef SLP
    return numFrames;

fail_alloc:
    /* All cleanup pointers are initialized before the first possible jump. */
    if (window) free(window);
    if (Bark2_freq) free(Bark2_freq);
    if (Bark2_bark) free(Bark2_bark);
    if (Barkno) free(Barkno);
    if (zb) free(zb);
    if (MinExcdB) free(MinExcdB);
    if (a0) free(a0);
    if (Hweight) free(Hweight);
    if (fft_buf) free(fft_buf);
    if (hBPi_all) free(hBPi_all);
    if (dataIn) free(dataIn);
    if (TempIn_re) free(TempIn_re);
    if (TempIn_im) free(TempIn_im);
    if (Slopes) free(Slopes);
    if (whichL) free(whichL);
    if (Lg) free(Lg);
    if (LdB) free(LdB);
    if (S2) free(S2);
    if (wzF) free(wzF);
    if (wzC) free(wzC);
    if (t_ifft_buf) { for (int t = 0; t < nThreads; t++) if (t_ifft_buf[t]) free(t_ifft_buf[t]); free(t_ifft_buf); }
    if (t_ei_k) { for (int t = 0; t < nThreads; t++) if (t_ei_k[t]) free(t_ei_k[t]); free(t_ei_k); }
    if (t_etmp_abs) { for (int t = 0; t < nThreads; t++) if (t_etmp_abs[t]) free(t_etmp_abs[t]); free(t_etmp_abs); }
    if (t_fft_etmp) { for (int t = 0; t < nThreads; t++) if (t_fft_etmp[t]) free(t_fft_etmp[t]); free(t_fft_etmp); }
    if (t_fft_plan) { for (int t = 0; t < nThreads; t++) if (t_fft_plan[t]) destroy_cfft_plan(t_fft_plan[t]); free(t_fft_plan); }
    destroy_cfft_plan(fft_plan);
    if (resampledAudio) free(resampledAudio);
    #undef HW
    #undef HBPI
    #undef SLP
    return computationError;
}
