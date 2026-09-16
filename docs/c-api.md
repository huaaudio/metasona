<!-- SPDX-FileCopyrightText: 2026 MetaSona contributors -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# C API: ABI version 1

Include the single public header:

```c
#include <metasona/metasona.h>
```

The header is C++-safe through `extern "C"`. Public functions are explicitly
exported; other native symbols are hidden. Installed CMake packages provide the
target `MetaSona::metasona`:

```cmake
find_package(MetaSona 0.1 REQUIRED)
target_link_libraries(my_program PRIVATE MetaSona::metasona)
```

When manually linking a static Windows build instead of using the exported
target, define `MS_STATIC` before including the header.

## ABI constants

| Constant | Value | Meaning |
|---|---:|---|
| `MS_ABI_VERSION` | `1` | Binary interface contract |
| `MS_SAMPLE_RATE_HZ` | `48000` | Only accepted native signal rate |
| `MS_THIRD_OCTAVE_BANDS` | `28` | Level-vector length |
| `MS_BARK_BANDS` | `240` | Specific-loudness length |

Call `ms_abi_version()` before using a dynamically discovered library and
reject a value other than the ABI your application was compiled against.
`ms_version_string()` returns a process-lifetime, read-only UTF-8 semantic
version string.

## Status values

Every calculation and size query returns a `ms_status`. No library function
prints, calls `exit`, changes global locale, or reports through `errno`.

| Status | Meaning |
|---|---|
| `MS_OK` | Success |
| `MS_ERROR_NULL_POINTER` | A required pointer is null |
| `MS_ERROR_INVALID_ARGUMENT` | Shape, count, range, or option is invalid |
| `MS_ERROR_INVALID_SAMPLE_RATE` | A signal rate is not exactly 48 kHz |
| `MS_ERROR_INVALID_SOUND_FIELD` | Field is neither free nor diffuse |
| `MS_ERROR_NONFINITE_INPUT` | Input contains NaN or infinity |
| `MS_ERROR_INPUT_TOO_SHORT` | No valid analysis can be produced |
| `MS_ERROR_OUTPUT_TOO_SMALL` | Caller output capacity is insufficient |
| `MS_ERROR_SIZE_OVERFLOW` | Checked size arithmetic would overflow |
| `MS_ERROR_ALLOCATION` | Temporary allocation failed |
| `MS_ERROR_NUMERICAL` | A finite valid result could not be produced |

`ms_status_string(status)` returns a process-lifetime, read-only English
description. It is diagnostic text, not a stable machine-readable identifier.

ABI v1 represents `ms_status` as `int32_t` and `ms_sound_field` as `uint32_t`;
the named choices remain enum-style integer constants. Public Windows calls
use an explicit `__cdecl` convention. The shared-library soname ABI is `1`,
independent of the pre-1.0 semantic project version.

## Read-only axes

```c
const double *ms_third_octave_centres_hz(void);
const double *ms_bark_axis(void);
```

The returned arrays have 28 and 240 elements respectively, remain valid for
the process lifetime, and must not be freed or modified. The Bark axis is 0.1,
0.2, …, 24.0 Bark.

## Stationary loudness from levels

```c
ms_status ms_loudness_from_levels(
    const double *levels_db,
    size_t level_count,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity);
```

`level_count` must equal 28. Values are finite dB SPL in the order returned by
`ms_third_octave_centres_hz()`. Values for the first 11 bands, 25 through
250 Hz, must not exceed 120 dB SPL; the implementation rejects higher values
instead of extrapolating its low-band model data. `loudness_sone` is required.
To omit specific loudness, pass `NULL, 0`; otherwise capacity must be at least
240. Input is read-only and never used as scratch storage.

## Stationary loudness from pressure

```c
ms_status ms_loudness_stationary(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *loudness_sone,
    double *specific_sone_per_bark,
    size_t specific_capacity);
```

Pressure is calibrated Pa and `sample_rate_hz` must equal 48,000. At least 256
samples are required. The optional specific-output convention matches the
level API. The C library never resamples.

## Time-varying loudness

Query first:

```c
ms_status ms_loudness_time_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);
```

For nonempty 48 kHz input, the count is `ceil(sample_count / 96)`. Calculation:

```c
ms_status ms_loudness_time(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *loudness_sone,
    size_t loudness_capacity,
    double *specific_sone_per_bark,
    size_t specific_capacity,
    size_t *frame_count_written);
```

Both output arrays are required. Specific output is frame-major and needs
`frame_count * 240` doubles. `frame_count_written` is required and receives the
required count on success or `MS_ERROR_OUTPUT_TOO_SMALL`, and zero on other
errors. Neither output array is changed on a failed calculation. Check multiplication before
allocating in caller code as well. The C result has no timestamps. Outputs sample the causal trajectory at input
indices 0, 96, 192, …; Python labels them at nominal centres 1, 3, 5, … ms.
These labels do not imply symmetric analysis windows.

## Roughness

```c
ms_status ms_roughness_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);

ms_status ms_roughness_dw(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    double *roughness_asper,
    size_t roughness_capacity,
    size_t *frame_count_written);
```

The complete-frame count is
`1 + floor((sample_count - 9600) / 4800)`; input shorter than 9,600 samples is
rejected. Query the count before allocating. On success, `frame_count_written`
equals the query result and output has that many asper values. On
`MS_ERROR_OUTPUT_TOO_SMALL`, it receives the required count; other failures
leave it unchanged. Failed calculations do not change the result buffer.

This ABI entry point uses the supplied reference's tabulated Daniel–Weber
model and spectral conventions. It remains experimental and must not
be presented as a normatively validated implementation of the full model.

## Tonality

```c
ms_status ms_tonality_frame_count(
    size_t sample_count,
    uint32_t sample_rate_hz,
    size_t *frame_count);

ms_status ms_tonality_aures(
    const double *pressure_pa,
    size_t sample_count,
    uint32_t sample_rate_hz,
    ms_sound_field field,
    double *tonality,
    size_t tonality_capacity,
    size_t *frame_count_written);
```

The complete-frame count is
`1 + floor((sample_count - 12000) / 6000)`; shorter input is rejected.
Output is nonnegative and may exceed one. On insufficient capacity,
`frame_count_written` receives the required count; other failures leave it
unchanged. Failed calculations do not change the result buffer.

Both roughness and tonality adapt the supplied psychohelperc kernels to
SQAT revision `e6228b789fc9`. Their
native sample counts must fit in `int` as well as an addressable double array;
larger counts return `MS_ERROR_SIZE_OVERFLOW`. The public function signatures
are unchanged. Tonality now uses a 250 ms window and 125 ms hop;
call the frame-count query when allocating outputs.

## Sharpness

```c
ms_status ms_sharpness_din(
    const double *specific_sone_per_bark,
    size_t specific_count,
    double *sharpness_acum);
```

Input must contain exactly 240 finite, nonnegative values in sone/Bark. A zero
pattern returns zero acum.

## Ownership, concurrency, and failure

- Callers retain ownership of every input and output buffer.
- Inputs are `const`, are never mutated, and may be reused after return.
- Returned metadata pointers are owned by the library and are read-only.
- All calculation scratch storage is per call; concurrent calls with distinct
  output buffers are supported.
- Every input range must be disjoint from every output range, and separate
  writable output buffers must be mutually disjoint.
- Output contents are unspecified after an error unless the individual
  function contract above explicitly states otherwise. Inspect `ms_status`
  before consuming results.

A complete allocation/query example is in [`../examples/metasona_example.c`](../examples/metasona_example.c).
