<!-- SPDX-FileCopyrightText: 2026 MetaSona contributors -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Algorithm, calibration, framing, and units

This document describes MetaSona 0.1.0 as implemented. It is not a
substitute for the standards or research publications that define the model
families.

## Common signal contract

Native signal functions accept a finite, nonempty mono sequence of calibrated
instantaneous acoustic pressure in pascals at exactly 48,000 samples/s. The
Python boundary accepts integer rates from 8 to 192 kHz and converts them once
with SciPy 1.15.3's polyphase FIR `resample_poly`. It uses a Kaiser 5.0 window,
constant-zero padding, and the exact reduced rational ratio, such as `160/147`
for 44.1 to 48 kHz. Reduced factors above 4,096 are rejected before filter
construction to prevent pathological resource use. No kernel performs hidden
resampling.

The sound-pressure level convention is

```text
L = 20 log10(p_rms / 20 µPa) dB SPL.
```

MetaSona does not infer microphone sensitivity, ADC full scale, gain, or
digital dBFS. Users must convert samples to Pa. Resampling acts on pressure and
does not intentionally renormalize amplitude, although FIR edge transients and
pass-band response can affect short signals.

`free` and `diffuse` select the free-field or diffuse-field correction used by
loudness-dependent calculations. They do not spatialize, mix, or otherwise
transform the input.

## Stationary loudness

The 28-level entry point consumes nominal one-third-octave levels at 25, 31.5,
40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000,
1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, and 12500 Hz.
Levels are dB SPL relative to 20 µPa.
The first 11 input bands, 25 through 250 Hz, are limited to 120 dB SPL because
the implemented low-band correction data do not define a supported
extrapolation above that level.

The Zwicker-model kernel applies low-frequency,
hearing-threshold, transmission, sound-field, critical-band, and upper-slope
relations, then integrates the specific-loudness pattern over Bark. Specific
loudness is returned at 0.1-Bark spacing from 0.1 to 24.0 Bark in sone/Bark;
the integrated value is in sone.

The stationary signal entry point passes pressure through 28 causal filters,
each containing three second-order sections from the supplied C reference.
It averages squared filter outputs over the entire record (legacy TimeSkip=0),
with no Hann window or mean subtraction, and calls the level kernel. Filters
start at zero, so record length and startup transients matter. The reference's
absolute 1e-12 Pa² energy floor is retained. At least 256 samples are required
after conversion to 48 kHz. No padding or native resampling is performed.

The target is ISO 532-1:2017, corrected English version dated 2017-11. This
release has not established normative conformance.

## Time-varying loudness

The time-varying path is stateful within one call but re-entrant across calls.
It uses the same tabulated third-octave filters, three frequency-dependent
energy lowpasses, and core values sampled at 2 kHz. The nonlinear stage has
two coupled decay states with 5/15/75 ms constants and 24 interpolation
substeps. The final 3.5/70 ms lowpasses also use 24 substeps and combine with
weights 0.47/0.53.

Outputs sample the causal trajectory at input indices 0, 96, 192, …, giving
`ceil(sample_count / 96)` results. There is no padded future interval. Prefix
results therefore remain unchanged when more input is appended. Both output
arrays are committed only after the complete calculation succeeds.

These are **causal results**, not values from symmetric 2 ms analysis windows.
The C ABI returns no time axis. For presentation consistency, Python labels
result `i` at the nominal centre `(i + 0.5) × 2 ms`; the first label is 1 ms.
That label must not be interpreted as an acausal centre-aligned estimate.
Filter state starts from zero at every call, so onset and short-record results
include startup behavior.

## Daniel–Weber-targeted roughness

Roughness uses complete 200 ms frames (9,600 samples), a 100 ms hop, a
symmetric Blackman window and a 9,600-point DFT. A complete last frame is
included; a trailing incomplete frame is ignored. Python labels centres at
100, 200, … ms. Each call owns its scratch arrays and commits results only
on success.

The source-informed stages use tabulated Bark interpolation, hearing
thresholds, outer/middle-ear transfer, 47 excitation channels, rectified
envelopes, H2/H5/H16/H21/H42 modulation responses and channel weighting from
the supplied `roughness_dw.c`. The five response groups and the reference's
502 Hz cutoff (500 Hz on this DFT grid) are retained. The H2 group ends at
355 Hz. No log-Gaussian surrogate or fitted calibration multiplier remains.

DFT bin k consistently means k × 5 Hz. Analysis covers 20–20,000 Hz; the
legacy port's extra one-based Bark offsets are removed. Outside the ±0.5 Bark
pass region, excitation is thresholded at the corresponding neighboring
channel. A second whole-channel threshold is not applied. Envelope RMS/mean
defines modulation depth, capped at one. Signed two-channel-spaced
correlations are preserved and the empirical channel weight multiplies the
squared depth/correlation term: `R = 0.25 Σ g(z) (m × correlations)²`.

Normalizing envelopes before the modulation transform avoids large
intermediate energies without changing depth or correlation. The normalized
1 kHz/70 Hz AM anchor is approximately 1.007067 asper at total 60 dB SPL.
This is a source-informed engineering model; published/normative validation
remains open.

## Aures-1985-targeted tonality

Tonality uses complete 80 ms frames (3,840 samples) with a 40 ms hop and a Hann
window. Python labels the first frame at 40 ms. A record shorter than 80 ms is
rejected and an incomplete trailing frame is ignored.

The independent implementation detects locally prominent spectral peaks,
interpolates their frequencies and levels, measures intrinsic bandwidth at the
half-power threshold (`peak - 3.0103 dB`), estimates noise and secondary-tone
masking, and combines frequency, bandwidth, level, and loudness weights. The
loudness ratio uses the requested free- or diffuse-field convention. The
result is a nonnegative dimensionless tonality value; MetaSona does not
label it as a standardized “tonality unit.” The published aggregation scale
factor `1.09` is retained without fixture-specific refitting. Output is not
capped at one: a multitone frame may exceed one because the result is neither
a probability nor a percentage.

Detection, masking, and tonal-energy removal consistently use the same
seven-line ±3-bin tonal group. Intrinsic bandwidth removes the independently
derived `1.44092`-bin full half-power width of the Hann analysis window. An
exact analytic symmetric-Hann response, expressed as three Dirichlet terms,
corrects parabolic frequency and peak-level bias without zero padding or a
copied correction table. Before estimating noise loudness, bounded sinusoidal
least-squares fits remove the detected narrowband tones from the frame; the
same seven-line cores remain excluded from the residual spectrum.

The interpreted Aures/Terhardt spectral domain is 20–5,000 Hz. With 12.5 Hz
bins and the implementation's required ±3-bin prominence guard, eligible peak
centres are 62.5–4,962.5 Hz. Tones at 5 kHz or above are therefore outside this
version's detector scope; the ±3-bin guard is an implementation choice and an
explicit interpretation uncertainty, not a claim about a normative cutoff.

Only combined tonality is public in ABI version 1. Component weights and
summary statistics are deliberately not exposed because they have not yet
completed independent cross-implementation validation.

## DIN 45692:2009-08-targeted sharpness weighting

Sharpness targets the weighting defined by
[DIN 45692:2009-08](https://www.dinmedia.de/en/standard/din-45692/117635111).
It consumes exactly one nonnegative 240-point specific-loudness pattern
at 0.1 through 24.0 Bark. It applies the DIN 45692 high-Bark weighting and
returns acum. A zero pattern returns zero. The API does not infer or recompute
specific loudness, accept a time series, or apply a sound-field correction.

Conformance has not been certified.

## Determinism and concurrency

Model tables are read-only. Per-call filter state, FFT buffers, and temporary
results are locally allocated, and the library has no mutable process-global
calculation state. Calls can be made concurrently when callers supply distinct
output buffers. The implementation disables floating-point contraction on
GCC/Clang builds and requests precise floating-point behavior from MSVC, but
last-bit differences across compilers and math libraries remain possible.
