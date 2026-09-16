// SPDX-License-Identifier: GPL-3.0-only
// Derived from SQAT: Sergio Aguirre and Gil Felix Greco, September 2026.
// MetaSona adaptation: Jiahua Zhang, September 2026.
// Port of il_find_narrowband in SQAT Tonality_Aures1985.m, revision
// e6228b789fc9a22251314f95678b9b1e08e60c55. See THIRD_PARTY.md and LICENSES/.
#ifndef METASONA_TONALITY_NARROWBAND_H
#define METASONA_TONALITY_NARROWBAND_H

/* SQAT implementation choices: at most 20 search attempts; a three-bin
 * power average; noise floors from the lower decile of half a critical band.
 * Width is the central 90% power span divided by 0.9, not a fitted gain. */
#define AURES_MAX_REGIONS 20
typedef struct aures_region {
    int first, last;
    double frequency, width, level, floor_left, floor_right;
} aures_region;

static double region_bridge(const aures_region *region, int bin)
{
    int span = region->last - region->first;
    if (span < 1) span = 1;
    return region->floor_left + (region->floor_right - region->floor_left)
        * (double)(bin - region->first) / (double)span;
}

static double region_quantile(const double *cdf, int first, int last, double df, double q)
{
    int k = first + 1;
    while (k < last && cdf[k] < q) ++k;
    return interp_xy(cdf[k-1], (double)(k-1)*df, cdf[k], (double)k*df, q);
}

/* scratch has seven count-sized double arrays; taken has count bytes.
 * They are allocated once per analysis call, reused for every frame. */
static int find_narrowband(const double *spl, int count, double df,
    double *scratch, unsigned char *taken, aures_region regions[AURES_MAX_REGIONS])
{
    double *average = scratch, *search = scratch + count;
    double *work = scratch + 2*count, *bark = scratch + 3*count;
    double *sorted = scratch + 4*count, *power = scratch + 5*count;
    double *cdf = scratch + 6*count;
    int found = 0;
    memset(taken, 0, (size_t)count);
    for (int k = 0; k < count; ++k) {
        average[k] = spl[k];
        if (k > 0 && k < count-1) {
            double sum = pow(10.0, spl[k-1]/10.0) + pow(10.0, spl[k]/10.0)
                       + pow(10.0, spl[k+1]/10.0);
            average[k] = 10.0*log10(sum) - 10.0*log10(3.0);
        }
        if (!isfinite(average[k])) return -1;
        search[k] = k*df < AURES_MIN_FREQUENCY || k*df > AURES_MAX_FREQUENCY
            ? -100.0 : average[k];
        work[k] = spl[k];
        bark[k] = fq_to_bark(k*df);
    }
    for (int attempt = 0; attempt < AURES_MAX_REGIONS; ++attempt) {
        int peak = 0, left, right, left_taken = 0, right_taken = 0;
        int begin, end, first, last;
        double maximum, fc, width, cbw, own = 0.0, low = 0.0, high = 0.0;
        double old_power = 0.0, floor_power = 0.0, total_power = 0.0;
        aures_region region;
        for (int k = 1; k < count; ++k) if (search[k] > search[peak]) peak = k;
        maximum = search[peak];
        if (maximum < 0.0) break;
        left = peak-1;
        for (;;) {
            if (left < 0) { left = 0; break; }
            if (average[left] + 3.0 < maximum) break;
            if (taken[left]) { left_taken = 1; break; }
            --left;
        }
        right = peak+1;
        for (;;) {
            if (right >= count) { right = count-1; break; }
            if (average[right] + 3.0 < maximum) break;
            if (taken[right]) { right_taken = 1; break; }
            ++right;
        }
        begin = left - round_to_int(0.5*(left*df < 500.0 ? 100.0 : 0.2*left*df)/df);
        if (begin < 0) begin = 0;
        first = left;
        if (!left_taken) {
            int n = left-begin+1, rank = (int)floor(0.1*n)-1;
            if (rank < 1) rank = 1; /* SQAT uses one-based order statistics. */
            memcpy(sorted, average+begin, (size_t)n*sizeof(double));
            qsort(sorted, (size_t)n, sizeof(double), compare_double_ascending);
            first = begin;
            for (int k = left; k >= begin; --k) {
                if (taken[k]) { first = k+1; break; }
                if (average[k] <= sorted[rank-1]) { first = k; break; }
            }
        }
        end = right + round_to_int(0.5*(right*df < 500.0 ? 100.0 : 0.2*right*df)/df);
        if (end >= count) end = count-1;
        last = right;
        if (!right_taken) {
            int n = end-right+1, rank = (int)floor(0.1*n)-1;
            if (rank < 1) rank = 1;
            memcpy(sorted, average+right, (size_t)n*sizeof(double));
            qsort(sorted, (size_t)n, sizeof(double), compare_double_ascending);
            last = end;
            for (int k = right; k <= end; ++k) {
                if (taken[k]) { last = k-1; break; }
                if (average[k] <= sorted[rank-1]) { last = k; break; }
            }
        }
        if (first < 0 || last >= count || first > last) return -1;
        for (int k = first; k <= last; ++k) { taken[k] = 1; search[k] = -100.0; }
        region.first = first; region.last = last;
        region.floor_left = work[first]; region.floor_right = work[last];
        fc = sqrt((left*df)*(right*df));
        width = (right-left)*df;
        cbw = 25.0 + 75.0*pow(1.0 + 1.4*pow(fc/1000.0, 2.0), 0.69);
        if (fc <= 0.0 || width <= 0.0 || width >= cbw) continue;
        {
            double z = fq_to_bark(fc);
            int low_count = 0, high_count = 0;
            for (int k = 0; k < count; ++k) {
                double p = pow(10.0, work[k]/10.0);
                if (bark[k] >= z-0.5 && bark[k] <= z+0.5) own += p;
                if (bark[k] >= z-1.5 && bark[k] < z-0.5) { low += p; ++low_count; }
                if (bark[k] > z+0.5 && bark[k] <= z+1.5) { high += p; ++high_count; }
            }
            if (!isfinite(own) || !isfinite(low) || !isfinite(high)) return -1;
            if (!low_count || !high_count) continue;
            if (10.0*log10(own)-10.0*log10(low) < AURES_TONE_THRESHOLD
                || 10.0*log10(own)-10.0*log10(high) < AURES_TONE_THRESHOLD) continue;
        }
        for (int k = left; k <= right; ++k) old_power += pow(10.0, work[k]/10.0);
        for (int k = first; k <= last; ++k) {
            double bridge = pow(10.0, region_bridge(&region, k)/10.0);
            if (k >= left && k <= right) floor_power += bridge;
            power[k] = fmax(pow(10.0, work[k]/10.0)-bridge, 0.0);
            total_power += power[k];
        }
        if (!isfinite(old_power) || !isfinite(floor_power) || !isfinite(total_power)) return -1;
        if (old_power <= floor_power || total_power <= 0.0 || last == first) continue;
        {
            double cumulative = 0.0;
            for (int k = first; k <= last; ++k) {
                cumulative += power[k];
                /* SQAT perturbation ensures strict monotonicity for interp1. */
                cdf[k] = cumulative/total_power + (k-first+1)*1e-12;
            }
        }
        region.frequency = fc;
        region.width = (region_quantile(cdf, first, last, df, 0.95)
                       -region_quantile(cdf, first, last, df, 0.05))/0.9;
        region.level = 10.0*log10(old_power-floor_power);
        if (!isfinite(region.width) || region.width <= 0.0 || !isfinite(region.level)) return -1;
        for (int k = first; k <= last; ++k) work[k] = region_bridge(&region, k);
        regions[found++] = region;
    }
    return found;
}
#endif
