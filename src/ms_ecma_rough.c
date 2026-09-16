// SPDX-License-Identifier: GPL-3.0-only
// MetaSona author: Jiahua Zhang, 2026.
// Adapted from SQAT/RefMap Roughness_ECMA418_2, Mike JB Lotinga and
// Matt Torjussen, University of Salford. See NOTICE.
// ECMA-418-2:2025 Section 7, without the optional entropy weighting.

#include "ms_ecma.h"
#include "third_party/pocketfft/pocketfft.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ROUGH_BLOCK 16384u
#define ROUGH_HOP 4096u
#define ROUGH_BINS 257u
#define ROUGH_DF (1500.0/512.0)

typedef struct { size_t bin; double prominence; } rough_peak;

static double *rough_alloc(size_t count)
{
    if (count > (size_t)PTRDIFF_MAX/sizeof(double)) return NULL;
    return (double *)calloc(count,sizeof(double));
}

static int compare_double(const void *a, const void *b)
{
    double x=*(const double *)a, y=*(const double *)b;
    return (x>y)-(x<y);
}

static int compare_peak(const void *a, const void *b)
{
    const rough_peak *x=(const rough_peak *)a, *y=(const rough_peak *)b;
    if (x->prominence != y->prominence)
        return (x->prominence<y->prominence)-(x->prominence>y->prominence);
    return (x->bin>y->bin)-(x->bin<y->bin);
}

static double refine_rate(const double *power, size_t k)
{
    static const double error[34]={0,.0457,.0907,.1346,.1765,.2157,.2515,
        .2828,.3084,.3269,.3364,.3348,.3188,.2844,.2259,.1351,0,
        -.1351,-.2259,-.2844,-.3188,-.3348,-.3364,-.3269,-.3084,
        -.2828,-.2515,-.2157,-.1765,-.1346,-.0907,-.0457,0,0};
    double beta[34], denominator=power[k-1u]-2.0*power[k]+power[k+1u];
    double estimate;
    size_t j, best=0u, upper;
    /* A flat-topped peak has no unique parabolic vertex. */
    if (denominator==0.0) return (double)k*ROUGH_DF;
    estimate=((double)k+.5*(power[k-1u]-power[k+1u])/denominator)*ROUGH_DF;
    for (j=0u;j<34u;++j) {
        beta[j]=(floor(estimate/ROUGH_DF)+(double)j/32.0)*ROUGH_DF-estimate-error[j];
        if (j<33u && fabs(beta[j])<fabs(beta[best])) best=j;
    }
    upper=best>0u && beta[best]*beta[best-1u]<0.0 ? best : best+1u;
    return estimate+error[upper-1u]-(error[upper]-error[upper-1u])*
        beta[upper-1u]/(beta[upper]-beta[upper-1u]);
}

static double rough_weight(double rate, double maximum, double a, double b)
{
    double x=a*(rate/maximum-maximum/rate);
    return pow(1.0+x*x,-b);
}

static double rough_estimate(const double *power, size_t band)
{
    rough_peak peaks[126];
    double rates[10], amps[10], max_peak=0.0, best_energy=-1.0;
    double f=ms_ecma_centre(band)/1000.0, lf=log2(f);
    double scale=1.0/(1.0+(f<1.0?.356:.8024)*pow(fabs(lf),f<1.0?.8049:.9333));
    double maximum=72.6937*(1.0-1.1739*exp(-5.4583*f));
    double hi_b=.2471+(lf>=-3.4253?.0129*(lf+3.4253)*(lf+3.4253):0.0);
    size_t k=3u, n=0u, i, j, used=0u, selected[10], best_set[10];
    size_t best_count=0u, best_fund=0u;
    while (k<255u) {
        size_t end=k;
        if (power[k]<=power[k-1u]) { ++k; continue; }
        while (end<255u && power[end+1u]==power[k]) ++end;
        if (end<255u && power[end]>power[end+1u]) {
            size_t pos=(k+end)/2u, p;
            double left=power[pos], right=power[pos];
            for (p=k;p>2u;) { --p; if (power[p]>power[pos]) break; left=fmin(left,power[p]); }
            for (p=end+1u;p<256u;++p) { if (power[p]>power[pos]) break; right=fmin(right,power[p]); }
            peaks[n].bin=pos; peaks[n].prominence=power[pos]-fmax(left,right); ++n;
        }
        k=end+1u;
    }
    if (n==0u) return 0.0;
    qsort(peaks,n,sizeof(peaks[0]),compare_peak);
    if (n>10u) n=10u;
    for (i=0u;i<n;++i) max_peak=fmax(max_peak,power[peaks[i].bin]);
    /* Restore ascending frequency order for deterministic harmonic ties. */
    for (i=1u;i<n;++i) {
        rough_peak v=peaks[i]; j=i;
        while (j>0u && peaks[j-1u].bin>v.bin) { peaks[j]=peaks[j-1u]; --j; }
        peaks[j]=v;
    }
    for (i=0u;i<n;++i) {
        k=peaks[i].bin;
        if (power[k]<=.05*max_peak) continue;
        rates[used]=refine_rate(power,k);
        amps[used]=(power[k-1u]+power[k]+power[k+1u])*scale;
        if (rates[used]<=ROUGH_DF) amps[used]=0.0;
        else if (rates[used]>maximum) amps[used]*=rough_weight(rates[used],maximum,1.2822,hi_b);
        ++used;
    }
    for (i=0u;i<used;++i) {
        size_t ratios[10], count=0u;
        double energy=0.0;
        for (j=0u;j<used;++j) ratios[j]=(size_t)floor(rates[j]/rates[i]+.5);
        for (j=0u;j<used;++j) {
            size_t p, closest=j;
            double error;
            if (ratios[j]==0u) continue;
            for (p=0u;p<j;++p) if (ratios[p]==ratios[j]) break;
            if (p<j) continue;
            error=fabs(rates[j]/((double)ratios[j]*rates[i])-1.0);
            for (p=j+1u;p<used;++p) if (ratios[p]==ratios[j]) {
                double candidate=fabs(rates[p]/((double)ratios[p]*rates[i])-1.0);
                if (candidate<error) { error=candidate; closest=p; }
            }
            if (error<.04) { selected[count++]=closest; energy+=amps[closest]; }
        }
        if (energy>best_energy) {
            best_energy=energy; best_count=count; best_fund=i;
            memcpy(best_set,selected,count*sizeof(size_t));
        }
    }
    if (best_count==0u || best_energy<=0.0) return 0.0;
    {
        size_t strongest=best_set[0];
        double centre=0.0, result, fundamental=rates[best_fund];
        for (i=0u;i<best_count;++i) {
            j=best_set[i]; centre+=rates[j]*amps[j];
            if (amps[j]>amps[strongest]) strongest=j;
        }
        centre/=best_energy+(double)best_count*1e-12;
        result=best_energy*(1.0+.1*pow(fabs(centre-rates[strongest]),.749));
        if (fundamental<=ROUGH_DF) return 0.0;
        if (fundamental<=maximum) result*=rough_weight(fundamental,maximum,.7066,1.0967-.064*lf);
        return result>=.074376 ? result : 0.0;
    }
}

static double pchip_endpoint(double h0, double h1, double d0, double d1)
{
    double slope=((2.0*h0+h1)*d0-h0*d1)/(h0+h1);
    if (slope*d0<=0.0) return 0.0;
    if (d0*d1<=0.0 && fabs(slope)>3.0*fabs(d0)) return 3.0*d0;
    return slope;
}

static void interpolate_rough(const double *est, const double *times, size_t frames,
    size_t nout, double *slopes, double *out)
{
    size_t b, i, j;
    for (b=0u;b<MS_ECMA_BANDS;++b) {
        double h0=times[1]-times[0], h1=times[2]-times[1];
        double d0=(est[MS_ECMA_BANDS+b]-est[b])/h0;
        double d1=(est[2u*MS_ECMA_BANDS+b]-est[MS_ECMA_BANDS+b])/h1;
        slopes[0]=pchip_endpoint(h0,h1,d0,d1);
        for (i=1u;i+1u<frames;++i) {
            double w1,w2;
            h0=times[i]-times[i-1u]; h1=times[i+1u]-times[i];
            d0=(est[i*MS_ECMA_BANDS+b]-est[(i-1u)*MS_ECMA_BANDS+b])/h0;
            d1=(est[(i+1u)*MS_ECMA_BANDS+b]-est[i*MS_ECMA_BANDS+b])/h1;
            w1=2.0*h1+h0; w2=h1+2.0*h0;
            slopes[i]=d0*d1<=0.0 ? 0.0 : (w1+w2)/(w1/d0+w2/d1);
        }
        slopes[frames-1u]=pchip_endpoint(h1,h0,d1,d0);
        j=0u;
        for (i=0u;i<nout;++i) {
            double t=(double)i/50.0, h, s, v, y0, y1;
            while (j+2u<frames && t>times[j+1u]) ++j;
            h=times[j+1u]-times[j]; s=(t-times[j])/h;
            y0=est[j*MS_ECMA_BANDS+b]; y1=est[(j+1u)*MS_ECMA_BANDS+b];
            v=(2.0*s*s*s-3.0*s*s+1.0)*y0+(s*s*s-2.0*s*s+s)*h*slopes[j]
              +(-2.0*s*s*s+3.0*s*s)*y1+(s*s*s-s*s)*h*slopes[j+1u];
            out[i*MS_ECMA_BANDS+b]=fmax(0.0,v);
        }
    }
}

ms_status ms_ecma_rough_count(size_t count, size_t *frames)
{
    if (frames==NULL) return MS_ERROR_NULL_POINTER;
    *frames=0u;
    if (count<15360u) return MS_ERROR_INPUT_TOO_SHORT;
    if (count>(size_t)PTRDIFF_MAX/sizeof(double)-ROUGH_BLOCK) return MS_ERROR_SIZE_OVERFLOW;
    *frames=count/960u+1u;
    return MS_OK;
}

ms_status MS_CALL ms_ecma_roughness_frame_count(size_t count, uint32_t rate, size_t *frames)
{
    if (frames==NULL) return MS_ERROR_NULL_POINTER;
    *frames=0u;
    if (rate!=MS_SAMPLE_RATE_HZ) return MS_ERROR_INVALID_SAMPLE_RATE;
    return ms_ecma_rough_count(count,frames);
}

ms_status MS_CALL ms_roughness_ecma(const double *input, size_t count,
    uint32_t rate, ms_sound_field field, double *out, size_t capacity, size_t *written)
{
    if (written==NULL) return MS_ERROR_NULL_POINTER;
    *written=0u;
    if (rate!=MS_SAMPLE_RATE_HZ) return MS_ERROR_INVALID_SAMPLE_RATE;
    if (!ms_is_valid_field(field)) return MS_ERROR_INVALID_SOUND_FIELD;
    return ms_ecma_rough(input,count,field,out,capacity,written);
}

ms_status ms_ecma_rough(const double *input, size_t count, ms_sound_field field,
    double *out, size_t capacity, size_t *written)
{
    double *ear=NULL, *bandpass=NULL, *powers=NULL, *loud=NULL, *work=NULL;
    double *est=NULL, *times=NULL, *slopes=NULL;
    double envelope[1024], window[512];
    double average[MS_ECMA_BANDS*ROUGH_BINS], summed[ROUGH_BINS];
    double clip[ROUGH_BINS], median[254];
    cfft_plan hilbert=NULL, spectrum=NULL;
    size_t padded, nout, frames, cells, b, f, k, start;
    ms_status status;
    if (written==NULL || out==NULL) return MS_ERROR_NULL_POINTER;
    status=ms_ecma_rough_count(count,&nout); *written=nout;
    if (status!=MS_OK) return status;
    if (capacity/MS_ECMA_BANDS<nout) return MS_ERROR_OUTPUT_TOO_SMALL;
    frames=count/ROUGH_HOP+1u+(count%ROUGH_HOP!=0u ? 1u : 0u);
    if (!ms_checked_multiply(frames,MS_ECMA_BANDS*ROUGH_BINS,&cells) ||
        cells>(size_t)PTRDIFF_MAX/sizeof(double)) return MS_ERROR_SIZE_OVERFLOW;
    status=ms_ecma_prepare(input,count,ROUGH_BLOCK,ROUGH_HOP,0,2025u,field,&ear,&padded);
    if (status!=MS_OK) return status;
    bandpass=rough_alloc(padded); powers=rough_alloc(cells);
    loud=rough_alloc(frames*MS_ECMA_BANDS); work=rough_alloc(2u*ROUGH_BLOCK);
    est=rough_alloc(frames*MS_ECMA_BANDS); times=rough_alloc(frames); slopes=rough_alloc(frames);
    hilbert=make_cfft_plan(ROUGH_BLOCK); spectrum=make_cfft_plan(512u);
    if (!bandpass || !powers || !loud || !work || !est || !times || !slopes || !hilbert || !spectrum) {
        status=MS_ERROR_ALLOCATION; goto done;
    }
    for (k=0u;k<512u;++k) window[k]=(.5-.5*cos(2.0*3.14159265358979323846*(double)k/512.0))/sqrt(.375);
    for (f=0u;f<frames;++f) times[f]=(double)(f+1u==frames?count:f*ROUGH_HOP)/48000.0;
    for (b=0u;b<MS_ECMA_BANDS;++b) {
        status=ms_ecma_band(ear,padded,b,bandpass);
        if (status!=MS_OK) goto done;
        for (f=0u;f<frames;++f) {
            double rms=0.0, energy=0.0, basis, scale;
            double *power=powers+(f*MS_ECMA_BANDS+b)*ROUGH_BINS;
            start=f+1u==frames?count:f*ROUGH_HOP;
            for (k=0u;k<ROUGH_BLOCK;++k) {
                double x=bandpass[start+k]; work[2u*k]=x; work[2u*k+1u]=0.0;
                if (x>0.0) rms+=x*x;
            }
            basis=ms_ecma_basis(sqrt(2.0*rms/(double)ROUGH_BLOCK),b,2025u);
            if (!ms_is_finite(basis)) { status=MS_ERROR_NUMERICAL; goto done; }
            loud[f*MS_ECMA_BANDS+b]=basis;
            if (basis==0.0) continue;
            if (cfft_forward(hilbert,work,1.0)) { status=MS_ERROR_NUMERICAL; goto done; }
            for (k=1u;k<ROUGH_BLOCK/2u;++k) { work[2u*k]*=2.0; work[2u*k+1u]*=2.0; }
            memset(work+ROUGH_BLOCK+2u,0,(ROUGH_BLOCK-2u)*sizeof(double));
            if (cfft_backward(hilbert,work,1.0/(double)ROUGH_BLOCK)) { status=MS_ERROR_NUMERICAL; goto done; }
            for (k=0u;k<512u;++k) {
                double x=hypot(work[64u*k],work[64u*k+1u])*window[k];
                envelope[2u*k]=x; envelope[2u*k+1u]=0.0; energy+=x*x;
            }
            if (cfft_forward(spectrum,envelope,1.0)) { status=MS_ERROR_NUMERICAL; goto done; }
            scale=energy>0.0?basis*basis/energy:0.0;
            for (k=0u;k<ROUGH_BINS;++k) {
                double re=envelope[2u*k], im=envelope[2u*k+1u];
                power[k]=scale*(re*re+im*im);
            }
        }
    }
    for (f=0u;f<frames;++f) {
        double max_loud=0.0, peak_clip=0.0, med;
        for (b=0u;b<MS_ECMA_BANDS;++b) max_loud=fmax(max_loud,loud[f*MS_ECMA_BANDS+b]);
        if (max_loud==0.0) continue;
        memset(summed,0,sizeof(summed));
        for (b=0u;b<MS_ECMA_BANDS;++b) for (k=0u;k<ROUGH_BINS;++k) {
            size_t ix=(f*MS_ECMA_BANDS+b)*ROUGH_BINS+k;
            double value=powers[ix];
            if (b>0u && b+1u<MS_ECMA_BANDS) value=(powers[ix-ROUGH_BINS]+value+powers[ix+ROUGH_BINS])/3.0;
            value/=max_loud;
            if (!ms_is_finite(value)) { status=MS_ERROR_NUMERICAL; goto done; }
            average[b*ROUGH_BINS+k]=value; summed[k]+=value;
        }
        memcpy(median,summed+2u,sizeof(median)); qsort(median,254u,sizeof(double),compare_double);
        med=.5*(median[126]+median[127])+1e-10;
        for (k=0u;k<ROUGH_BINS;++k) {
            clip[k]=.0856*summed[k]/med*fmin(.1891*exp(.012*(double)k),1.0);
            if (k>=2u && k<256u) peak_clip=fmax(peak_clip,clip[k]);
        }
        for (k=0u;k<ROUGH_BINS;++k) clip[k]=clip[k]>=.05*peak_clip?fmin(fmax(clip[k]-.1407,0.0),1.0):0.0;
        for (b=0u;b<MS_ECMA_BANDS;++b) {
            for (k=0u;k<ROUGH_BINS;++k) average[b*ROUGH_BINS+k]*=clip[k];
            est[f*MS_ECMA_BANDS+b]=rough_estimate(average+b*ROUGH_BINS,b);
            if (!ms_is_finite(est[f*MS_ECMA_BANDS+b])) { status=MS_ERROR_NUMERICAL; goto done; }
        }
    }
    interpolate_rough(est,times,frames,nout,slopes,out);
    for (f=0u;f<nout;++f) {
        double sum=0.0, square=0.0, ratio, exponent;
        for (b=0u;b<MS_ECMA_BANDS;++b) {
            double x=out[f*MS_ECMA_BANDS+b]; sum+=x; square+=x*x;
        }
        ratio=sum>0.0?sqrt(square*(double)MS_ECMA_BANDS)/sum:0.0;
        exponent=.37106*(tanh(1.6407*(ratio-2.5804))+1.0)*.5+.58449;
        for (b=0u;b<MS_ECMA_BANDS;++b) {
            size_t ix=f*MS_ECMA_BANDS+b;
            double x=.0180685*pow(out[ix],exponent);
            if (f>0u) {
                double previous=out[ix-MS_ECMA_BANDS];
                double d=exp(-1.0/(50.0*(x>=previous?.0625:.5)));
                x=(1.0-d)*x+d*previous;
            }
            if (!ms_is_finite(x)) { status=MS_ERROR_NUMERICAL; goto done; }
            out[ix]=x;
        }
    }
    status=MS_OK;
done:
    free(ear); free(bandpass); free(powers); free(loud); free(work);
    free(est); free(times); free(slopes);
    if (hilbert) destroy_cfft_plan(hilbert);
    if (spectrum) destroy_cfft_plan(spectrum);
    return status;
}
