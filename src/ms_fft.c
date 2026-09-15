// SPDX-License-Identifier: Apache-2.0

#include "ms_internal.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

int ms_next_power_of_two(size_t value, size_t *result)
{
    size_t power = 1u;
    if (result == NULL || value == 0u) {
        return 0;
    }
    while (power < value) {
        if (power > SIZE_MAX / 2u) {
            return 0;
        }
        power *= 2u;
    }
    *result = power;
    return 1;
}

static ms_complex ms_complex_multiply(ms_complex left, ms_complex right)
{
    ms_complex result;
    result.re = left.re * right.re - left.im * right.im;
    result.im = left.re * right.im + left.im * right.re;
    return result;
}

static ms_status ms_fft_power_two(ms_complex *values, size_t count, int inverse)
{
    size_t index;
    size_t reverse = 0u;
    size_t length;
    if (values == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (count < 2u || (count & (count - 1u)) != 0u) {
        return MS_ERROR_INVALID_ARGUMENT;
    }

    for (index = 1u; index < count; ++index) {
        size_t bit = count >> 1u;
        ms_complex temporary;
        while ((reverse & bit) != 0u) {
            reverse ^= bit;
            bit >>= 1u;
        }
        reverse ^= bit;
        if (index < reverse) {
            temporary = values[index];
            values[index] = values[reverse];
            values[reverse] = temporary;
        }
    }

    for (length = 2u; length <= count; length <<= 1u) {
        const double angle = (inverse ? 2.0 : -2.0) * MS_PI / (double)length;
        ms_complex step;
        size_t offset;
        step.re = cos(angle);
        step.im = sin(angle);
        for (offset = 0u; offset < count; offset += length) {
            ms_complex twiddle = {1.0, 0.0};
            const size_t half = length >> 1u;
            size_t element;
            for (element = 0u; element < half; ++element) {
                const ms_complex even = values[offset + element];
                const ms_complex odd = ms_complex_multiply(
                    twiddle, values[offset + element + half]);
                values[offset + element].re = even.re + odd.re;
                values[offset + element].im = even.im + odd.im;
                values[offset + element + half].re = even.re - odd.re;
                values[offset + element + half].im = even.im - odd.im;
                twiddle = ms_complex_multiply(twiddle, step);
            }
        }
        if (length == count) {
            break;
        }
    }
    if (inverse) {
        for (index = 0u; index < count; ++index) {
            values[index].re /= (double)count;
            values[index].im /= (double)count;
        }
    }
    return MS_OK;
}

static ms_status ms_fft_bluestein(ms_complex *values, size_t count, int inverse)
{
    ms_complex *left;
    ms_complex *right;
    size_t convolution_size;
    size_t allocation_size;
    size_t index;
    ms_status status;
    if (count > (SIZE_MAX / 2u) + 1u
        || !ms_next_power_of_two(2u * count - 1u, &convolution_size)
        || !ms_checked_multiply(convolution_size, sizeof(*left), &allocation_size)) {
        return MS_ERROR_SIZE_OVERFLOW;
    }
    left = (ms_complex *)calloc(1u, allocation_size);
    right = (ms_complex *)calloc(1u, allocation_size);
    if (left == NULL || right == NULL) {
        free(left);
        free(right);
        return MS_ERROR_ALLOCATION;
    }
    for (index = 0u; index < count; ++index) {
        const double position = (double)index;
        const double angle = MS_PI * position * position / (double)count;
        const double sign = inverse ? 1.0 : -1.0;
        ms_complex input_chirp = {cos(angle), sign * sin(angle)};
        ms_complex convolution_chirp = {cos(angle), -sign * sin(angle)};
        left[index] = ms_complex_multiply(values[index], input_chirp);
        right[index] = convolution_chirp;
        if (index != 0u) {
            right[convolution_size - index] = convolution_chirp;
        }
    }
    status = ms_fft_power_two(left, convolution_size, 0);
    if (status == MS_OK) {
        status = ms_fft_power_two(right, convolution_size, 0);
    }
    if (status == MS_OK) {
        for (index = 0u; index < convolution_size; ++index) {
            left[index] = ms_complex_multiply(left[index], right[index]);
        }
        status = ms_fft_power_two(left, convolution_size, 1);
    }
    if (status == MS_OK) {
        for (index = 0u; index < count; ++index) {
            const double position = (double)index;
            const double angle = MS_PI * position * position / (double)count;
            const double sign = inverse ? 1.0 : -1.0;
            ms_complex output_chirp = {cos(angle), sign * sin(angle)};
            values[index] = ms_complex_multiply(left[index], output_chirp);
            if (inverse) {
                values[index].re /= (double)count;
                values[index].im /= (double)count;
            }
        }
    }
    free(left);
    free(right);
    return status;
}

static ms_status ms_fft_transform(ms_complex *values, size_t count, int inverse)
{
    ms_status status;
    size_t index;
    if (values == NULL) {
        return MS_ERROR_NULL_POINTER;
    }
    if (count < 2u) {
        return MS_ERROR_INVALID_ARGUMENT;
    }
    for (index = 0u; index < count; ++index) {
        if (!ms_is_finite(values[index].re) || !ms_is_finite(values[index].im)) {
            return MS_ERROR_NUMERICAL;
        }
    }
    if ((count & (count - 1u)) == 0u) {
        status = ms_fft_power_two(values, count, inverse);
    } else {
        status = ms_fft_bluestein(values, count, inverse);
    }
    if (status != MS_OK) {
        return status;
    }
    for (index = 0u; index < count; ++index) {
        if (!ms_is_finite(values[index].re) || !ms_is_finite(values[index].im)) {
            return MS_ERROR_NUMERICAL;
        }
    }
    return MS_OK;
}

ms_status ms_fft_forward(ms_complex *values, size_t count)
{
    return ms_fft_transform(values, count, 0);
}

ms_status ms_fft_inverse(ms_complex *values, size_t count)
{
    return ms_fft_transform(values, count, 1);
}
