// SPDX-License-Identifier: GPL-3.0-only
// MetaSona author: Jiahua Zhang, 2026.
// Adapted from SQAT/RefMap Tonality_ECMA418_2 and Loudness_ECMA418_2,
// Mike JB Lotinga and Matt Torjussen, University of Salford. See NOTICE.
// ECMA-418-2:2025 Sections 6 and 8. Call-local buffers and FFT plans.

#include "ms_ecma.h"
#include "third_party/pocketfft/pocketfft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static double *new_doubles(size_t count)
{
    if (count > (size_t)PTRDIFF_MAX / sizeof(double)) return NULL;
    return (double *)calloc(count, sizeof(double));
}

ms_status MS_CALL ms_ecma_tonal_frame_count(size_t count, uint32_t rate, size_t *frames)
{
    if (frames==NULL) return MS_ERROR_NULL_POINTER;
    *frames=0u;
    if (rate!=MS_SAMPLE_RATE_HZ) return MS_ERROR_INVALID_SAMPLE_RATE;
    return ms_ecma_tonal_count(count,frames);
}

ms_status MS_CALL ms_ecma_tonal_analysis(const double *input, size_t count,
    uint32_t rate, ms_sound_field field, double *out, size_t capacity, size_t *written)
{
    if (written==NULL) return MS_ERROR_NULL_POINTER;
    *written=0u;
    if (rate!=MS_SAMPLE_RATE_HZ) return MS_ERROR_INVALID_SAMPLE_RATE;
    if (!ms_is_valid_field(field)) return MS_ERROR_INVALID_SOUND_FIELD;
    return ms_ecma_tonal(input,count,field,out,capacity,written);
}

static void smooth(double *values, size_t count)
{
    double d = exp(-1.0 / (187.5 * (6.0/224.0)));
    double gain = pow(1.0-d, 3.0)/(d+d*d);
    double a1 = -3.0*d, a2 = 3.0*d*d, a3 = -d*d*d;
    double z0 = 0.0, z1 = 0.0, z2 = 0.0;
    size_t i;
    for (i = 0u; i < count; ++i) {
        double x = values[i], y = z0;
        z0 = gain*d*x - a1*y + z1;
        z1 = gain*d*d*x - a2*y + z2;
        z2 = -a3*y;
        values[i] = y;
    }
}

static ms_status scaled_acf(const double *signal, size_t padded, size_t band,
    size_t block, size_t frames, size_t lags, rfft_plan plan, double *work,
    double *prefix, double *suffix, double *out)
{
    size_t f, i;
    size_t hop = block/4u, start = 8192u-block;
    (void)padded;
    for (f = 0u; f < frames; ++f) {
        double basis, fft_energy;
        memset(work, 0, 2u*block*sizeof(double));
        prefix[0] = 0.0;
        for (i = 0u; i < block; ++i) {
            double x = signal[start + f*hop + i];
            if (x < 0.0) x = 0.0;
            work[i] = x;
            prefix[i+1u] = prefix[i] + x*x;
        }
        basis = ms_ecma_basis(sqrt(2.0*prefix[block]/(double)block), band, 2025u);
        if (!ms_is_finite(basis)) return MS_ERROR_NUMERICAL;
        if (basis == 0.0) continue; /* output is zero-initialized */
        suffix[block] = 0.0;
        for (i = block; i > 0u; --i) {
            double x = work[i-1u];
            suffix[i-1u] = suffix[i] + x*x;
        }
        if (rfft_forward(plan, work, 1.0)) return MS_ERROR_NUMERICAL;
        /* PocketFFT real packing: DC, Re(1), Im(1), ..., Nyquist. */
        work[0] *= work[0];
        fft_energy=work[0];
        for (i = 1u; i < block; ++i) {
            double re = work[2u*i-1u], im = work[2u*i];
            work[2u*i-1u] = re*re + im*im;
            work[2u*i] = 0.0;
            fft_energy+=work[2u*i-1u];
        }
        work[2u*block-1u] *= work[2u*block-1u];
        fft_energy+=work[2u*block-1u];
        if (!ms_is_finite(fft_energy)) return MS_ERROR_NUMERICAL;
        if (rfft_backward(plan, work, 1.0/(double)(2u*block))) return MS_ERROR_NUMERICAL;
        for (i = 0u; i < lags; ++i) {
            if (i < 3u*block/4u && suffix[i]>0.0 && prefix[block-i]>0.0) {
                /* A zero-energy overlap has exactly zero autocorrelation.
                 * Otherwise FFT roundoff divided by epsilon can create a
                 * spurious component during the leading-zero transient.
                 * Rectification and Cauchy-Schwarz bound the normalized
                 * autocorrelation to [0,1]. */
                double normalized=work[i]/(sqrt(suffix[i])*sqrt(prefix[block-i])+1e-12);
                out[f*lags+i] = basis*fmin(1.0,fmax(0.0,normalized));
            }
        }
    }
    return MS_OK;
}

ms_status ms_ecma_tonal_count(size_t count, size_t *frames)
{
    size_t n;
    if (frames == NULL) return MS_ERROR_NULL_POINTER;
    *frames = 0u;
    if (count < 14592u) return MS_ERROR_INPUT_TOO_SHORT; /* 304 ms */
    if (count > (size_t)PTRDIFF_MAX/sizeof(double) - 32768u) return MS_ERROR_SIZE_OVERFLOW;
    n = (count+255u)/256u + 1u;
    if (n > (size_t)PTRDIFF_MAX / (5u*MS_ECMA_BANDS*sizeof(double)))
        return MS_ERROR_SIZE_OVERFLOW;
    *frames = n;
    return MS_OK;
}

ms_status ms_ecma_tonal(const double *input, size_t count, ms_sound_field field,
    double *out, size_t capacity, size_t *written)
{
    static const size_t low[4] = {0u,3u,16u,25u}, high[4] = {2u,15u,24u,52u};
    static const size_t blocks[4] = {8192u,4096u,2048u,1024u};
    static const double cs[4] = {18.21,12.14,417.54,962.68};
    static const double ds[4] = {.36,.36,.71,.69};
    double *ear = NULL, *bandpass = NULL, *acf = NULL, *average = NULL;
    double *work = NULL, *prefix = NULL, *suffix = NULL;
    double *raw = NULL, *interpolated = NULL;
    rfft_plan plan = NULL, spectral = NULL;
    size_t padded, nout, plane, group, band, i, f;
    ms_status status;
    if (written == NULL || out == NULL) return MS_ERROR_NULL_POINTER;
    status = ms_ecma_tonal_count(count, &nout);
    *written = nout;
    if (status != MS_OK) return status;
    plane = nout*MS_ECMA_BANDS;
    if (capacity < 5u*plane) return MS_ERROR_OUTPUT_TOO_SMALL;
    status = ms_ecma_prepare(input,count,8192u,2048u,1,2025u,field,&ear,&padded);
    if (status != MS_OK) return status;
    bandpass = new_doubles(padded);
    work = new_doubles(16384u);
    prefix = new_doubles(8193u);
    suffix = new_doubles(8193u);
    interpolated = new_doubles(4u*nout);
    spectral = make_rfft_plan(16384u);
    if (!bandpass || !work || !prefix || !suffix || !interpolated || !spectral) {
        status = MS_ERROR_ALLOCATION; goto done;
    }
    for (group = 0u; group < 4u; ++group) {
        size_t block = blocks[group], hop = block/4u, factor = block/1024u;
        size_t frames = (padded-8192u)/hop+1u;
        size_t left = group == 0u ? 0u : (group == 1u ? 1u : (group == 2u ? 15u : 25u));
        size_t right = group == 0u ? 4u : (group == 1u ? 17u : (group == 2u ? 25u : 52u));
        size_t neighbors = right-left+1u;
        double delta = hypot(81.9289,.1618*ms_ecma_centre(low[group]));
        double tau_start = fmax(.5/delta,.002);
        size_t lags = (size_t)floor(fmax(4.0/delta,tau_start+.001)*48000.0);
        size_t cells;
        if (!ms_checked_multiply(frames,lags,&cells) ||
            cells > (size_t)PTRDIFF_MAX/(neighbors*sizeof(double))) {
            status = MS_ERROR_SIZE_OVERFLOW; goto done;
        }
        acf = new_doubles(neighbors*cells);
        average = new_doubles(cells);
        raw = new_doubles(3u*frames);
        plan = make_rfft_plan(2u*block);
        if (!acf || !average || !raw || !plan) { status = MS_ERROR_ALLOCATION; goto done; }
        for (band = left; band <= right; ++band) {
            status = ms_ecma_band(ear,padded,band,bandpass);
            if (status != MS_OK) goto done;
            status = scaled_acf(bandpass,padded,band,block,frames,lags,plan,
                                work,prefix,suffix,acf+(band-left)*cells);
            if (status != MS_OK) goto done;
        }
        for (band = low[group]; band <= high[group]; ++band) {
            size_t before = band == 0u ? 0u : (band == 1u ? 1u : (band < 16u ? 2u : (band < 25u ? 1u : 0u)));
            size_t after = band < 2u ? 1u : before;
            size_t nband = before+after+1u, b;
            double frequency = ms_ecma_centre(band);
            double width = hypot(81.9289,.1618*frequency);
            double begin = fmax(.5/width,.002);
            size_t first = (size_t)ceil(begin*48000.0)-1u;
            size_t last = (size_t)floor(fmax(4.0/width,begin+.001)*48000.0)-1u;
            double win = (double)(last-first+1u);
            double g = cs[group]/pow(frequency,ds[group]);
            memset(average,0,cells*sizeof(double));
            for (b = band-before; b <= band+after; ++b)
                for (i = 0u; i < cells; ++i) average[i] += acf[(b-left)*cells+i]/(double)nband;
            for (f = 0u; f < frames; ++f) {
                double mean = 0.0, maximum = 0.0, spectral_energy = 0.0;
                size_t best = 0u;
                int temporal = band < 16u && f > 0u && f+1u < frames;
                double loudness = average[f*lags];
                if (temporal) loudness = (loudness+average[(f-1u)*lags]+average[(f+1u)*lags])/3.0;
                if (loudness==0.0) {
                    raw[f]=0.0; raw[frames+f]=0.0; raw[2u*frames+f]=0.0;
                    continue;
                }
                memset(work,0,16384u*sizeof(double));
                for (i = first; i <= last; ++i) {
                    double a = average[f*lags+i];
                    if (temporal) a = (a+average[(f-1u)*lags+i]+average[(f+1u)*lags+i])/3.0;
                    work[i] = a;
                    mean += a/win;
                }
                for (i = first; i <= last; ++i) work[i] -= mean;
                if (rfft_forward(spectral,work,1.0)) { status=MS_ERROR_NUMERICAL; goto done; }
                /* Real ACF has conjugate symmetry; use positive frequencies. */
                for (i = 0u; i <= 8192u; ++i) {
                    double power = i==0u ? work[0]*work[0] : (i==8192u ?
                        work[16383u]*work[16383u] :
                        work[2u*i-1u]*work[2u*i-1u]+work[2u*i]*work[2u*i]);
                    spectral_energy += power;
                    if (power > maximum) { maximum=power; best=i; }
                }
                if (!ms_is_finite(spectral_energy)) { status=MS_ERROR_NUMERICAL; goto done; }
                raw[f] = fmin(loudness,4.0*sqrt(maximum)/win);
                raw[frames+f] = loudness;
                raw[2u*frames+f] = (double)best*48000.0/16384.0;
            }
            for (i = 0u; i < nout; ++i) {
                size_t j = i/factor, next = j+1u < frames ? j+1u : j;
                double frac = (double)(i%factor)/(double)factor;
                double tone = raw[j]*(1.0-frac)+raw[next]*frac;
                double total = raw[frames+j]*(1.0-frac)+raw[frames+next]*frac;
                interpolated[i] = tone;
                interpolated[nout+i] = total;
                interpolated[2u*nout+i] = tone/(fmax(0.0,total-tone)+1e-12);
                out[2u*plane+i*MS_ECMA_BANDS+band] = raw[2u*frames+j]*(1.0-frac)+raw[2u*frames+next]*frac;
            }
            smooth(interpolated,nout);
            smooth(interpolated+nout,nout);
            smooth(interpolated+2u*nout,nout);
            for (i = 0u; i < nout; ++i) {
                double ratio = interpolated[2u*nout+i]/g;
                double reduction = ratio > .07 ? -expm1(-20.0*(ratio-.07)) : 0.0;
                double tone = fmax(0.0,reduction*interpolated[i]);
                out[i*MS_ECMA_BANDS+band] = tone;
                out[plane+i*MS_ECMA_BANDS+band] = fmax(0.0,interpolated[nout+i]-tone);
            }
        }
        free(acf); acf=NULL; free(average); average=NULL; free(raw); raw=NULL;
        destroy_rfft_plan(plan); plan=NULL;
    }
    for (i = 0u; i < nout; ++i) {
        double max_tone=0.0, noise_sum=0.0, max_total=0.0, ratio, q, exponent;
        for (band = 0u; band < MS_ECMA_BANDS; ++band) {
            double t=out[i*MS_ECMA_BANDS+band], n=out[plane+i*MS_ECMA_BANDS+band];
            max_tone=fmax(max_tone,t); noise_sum+=n; max_total=fmax(max_total,t+n);
        }
        ratio=max_tone/(noise_sum+1e-12);
        q=ratio > .003 ? -expm1(-35.0*(ratio-.003)) : 0.0;
        exponent=.2918/(max_total+1e-12)+.5459;
        for (band=0u;band<MS_ECMA_BANDS;++band) {
            size_t index=i*MS_ECMA_BANDS+band;
            double t=out[index], n=.5331*out[plane+index], scale=fmax(t,n);
            out[3u*plane+index]=(2.8758615/.9999043734252)*q*t;
            out[4u*plane+index]=scale > 0.0 ? scale*pow(pow(t/scale,exponent)+pow(n/scale,exponent),1.0/exponent) : 0.0;
        }
    }
    for (i=0u;i<5u*plane;++i) if (!ms_is_finite(out[i])) { status=MS_ERROR_NUMERICAL; goto done; }
    status=MS_OK;
done:
    free(ear); free(bandpass); free(acf); free(average); free(work);
    free(prefix); free(suffix); free(raw); free(interpolated);
    if (plan) destroy_rfft_plan(plan);
    if (spectral) destroy_rfft_plan(spectral);
    return status;
}
