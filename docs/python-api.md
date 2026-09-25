<!-- SPDX-FileCopyrightText: 2026 MetaSona contributors -->
<!-- SPDX-License-Identifier: Apache-2.0 -->

# Python API

All public names below are exported by `metasona`. Import is side-effect
free: the shared library is located, loaded, and checked only on the first
native numerical call.

## Types

`SoundField` is a string enum with `SoundField.FREE` (`"free"`) and
`SoundField.DIFFUSE` (`"diffuse"`). Public result dataclasses are frozen; every
array is a C-contiguous, read-only `numpy.float64` defensive copy.

| Result type | Fields |
|---|---|
| `LoudnessResult` | `loudness_sone`, `specific_loudness_sone_per_bark` `(240,)`, `bark_axis` `(240,)` |
| `TimeVaryingLoudnessResult` | `time_s` `(F,)`, `loudness_sone` `(F,)`, `specific_loudness_sone_per_bark` `(F,240)`, `bark_axis` `(240,)` |
| `RoughnessResult` | `time_s` `(F,)`, `roughness_asper` `(F,)` |
| `TonalityResult` | `time_s` `(F,)`, `tonality` `(F,)` |
| `SharpnessResult` | `sharpness_acum` |

## Calculations

```python
loudness_from_levels(
    levels_db_spl,
    *,
    sound_field=SoundField.FREE,
) -> LoudnessResult
```

`levels_db_spl` must have shape `(28,)`, be finite, and follow the 25 Hz to
12.5 kHz nominal-centre order returned by `third_octave_centres_hz()`. The
first 11 values, 25 through 250 Hz, must not exceed 120 dB SPL; larger values
would extrapolate beyond the implemented low-band model data and are rejected.

```python
stationary_loudness(
    pressure_pa,
    sample_rate_hz,
    *,
    sound_field=SoundField.FREE,
) -> LoudnessResult
```

The resampled signal must contain at least 256 samples. The complete input is
the stationary analysis interval.

```python
time_varying_loudness(
    pressure_pa,
    sample_rate_hz,
    *,
    sound_field=SoundField.FREE,
    time_skip_s=0.0,
) -> TimeVaryingLoudnessResult
```

Returns causal 2 ms results labelled at nominal frame centres: 0.001, 0.003,
0.005 seconds, and so on. These are presentation labels for causal samples
at 0, 2, 4, … ms. No padding is applied; the count is `ceil(samples / 96)`.

```python
roughness_daniel_weber(
    pressure_pa,
    sample_rate_hz,
    *,
    time_skip_s=0.0,
) -> RoughnessResult
```

Returns complete 200 ms frames with 100 ms hop, labelled 0.1, 0.2, … seconds.
This function is experimental. It uses the supplied reference's tabulated
model and spectral conventions; normative validation remains open.

```python
tonality_aures(
    pressure_pa,
    sample_rate_hz,
    *,
    sound_field=SoundField.FREE,
    time_skip_s=0.0,
) -> TonalityResult
```

Returns complete 250 ms frames with 125 ms hop, labelled 0.125, 0.250, … seconds.
The minimum input duration is 250 ms after resampling. Values are nonnegative
and may exceed one for multitone signals. Both metrics adapt the supplied
psychohelperc kernels to SQAT revision `e6228b789fc9`, and include the final
complete frame. See [implementation differences](../THIRD_PARTY.md).

```python
sharpness_din45692(specific_loudness_sone_per_bark) -> SharpnessResult
```

The input must be one nonnegative, finite pattern with shape `(240,)`.

## Metadata

- `native_version() -> str` loads the library and reports its semantic version.
- `third_octave_centres_hz() -> ndarray` returns the 28-value read-only nominal
  centre-frequency vector.
- `__version__` reports the Python package version.
- `NativeStatus` mirrors the stable ABI status codes.

## Validation and conversion

Signal inputs are converted to `float64`, then rejected unless they are
one-dimensional, nonempty, and finite. This conversion is an independent copy;
the caller's array is not mutated. Stereo or other multidimensional input is
never flattened implicitly.

`sample_rate_hz` must be an integer in `[8000, 192000]`; booleans, floats, and
numeric strings are rejected. Rates whose reduced ratio has a factor above
4,096 are rejected before filter allocation. Signals not already at 48 kHz are
converted with SciPy `scipy.signal.resample_poly`, an explicit Kaiser
5.0 window, constant-zero padding, and reduced integer factors. For example,
44.1 kHz uses `up=160`, `down=147`.

`time_skip_s` must be finite, nonnegative, shorter than the original input
duration, and leave at least one labelled output frame. Frames are retained
when `time_s >= time_skip_s`; the kernels still process the complete signal, so
skip does not reset their state.

## Rolling analysis

`RollingAnalyzer` adapts the batch signal metrics to applications that receive
calibrated mono pressure incrementally. It is a synchronous rolling-window
wrapper, not an audio-capture API or a stateful native streaming kernel.

```python
analyzer = RollingAnalyzer(
  sample_rate_hz=48_000,
  metrics=[
    RollingMetric.STATIONARY_LOUDNESS,
    RollingMetric.TONALITY_AURES,
  ],
  window_s=1.0,
  hop_s=0.2,
  sound_field=SoundField.FREE,
)

snapshots = analyzer.push(pressure_chunk_pa)
```

`metrics` is a nonempty iterable of unique `RollingMetric` values or their
exact strings:

| Metric | Result type |
|---|---|
| `STATIONARY_LOUDNESS` | `LoudnessResult` |
| `TIME_VARYING_LOUDNESS` | `TimeVaryingLoudnessResult` |
| `ROUGHNESS_DANIEL_WEBER` | `RoughnessResult` |
| `TONALITY_AURES` | `TonalityResult` |
| `SHARPNESS_DIN45692` | `SharpnessResult` |
| `LOUDNESS_ECMA` | `EcmaLoudnessResult` |
| `TONALITY_ECMA` | `EcmaTonalityResult` |
| `ROUGHNESS_ECMA` | `EcmaRoughnessResult` |

The analyzer rounds `window_s` and `hop_s` to input sample counts once. The
first snapshot ends after one complete window; later snapshots are separated
by one hop. `push` returns an empty tuple before that point and may return
multiple `RollingSnapshot` objects when a large chunk crosses several
boundaries. Snapshot start/end times are absolute on the stream timeline;
timestamps inside nested metric results remain relative to that window.

Chunks must be finite, nonempty, one-dimensional real arrays at the sample
rate fixed in the constructor. Input is copied. Completed windows are emitted
only without padding; there is no partial-window flush. `reset()` discards the
buffer and restarts timestamps at zero after a stream discontinuity.

Each snapshot contains an immutable `results` mapping. Sharpness reuses a
stationary-loudness calculation, and selected ECMA loudness and tonality share
one combined analysis. All due work runs synchronously and no snapshot is
dropped. Instances are not thread-safe; applications should serialize calls
and run expensive selections outside a real-time capture callback. As a
reference on a Ryzen 9 5950X, one-second windows took approximately 10 ms for
stationary loudness, 41 ms for Aures tonality, 138-158 ms for either roughness
model, and 615 ms for combined ECMA loudness/tonality. Choose a hop suitable
for the selected metrics and target hardware.

See the [complete rolling example](../examples/metasona_rolling.py).

## Exceptions and native loading

- `MetaSonaValidationError` reports an invalid Python argument before a
  native call.
- `NativeLibraryError` reports missing libraries, loading failures, an ABI
  mismatch, or a native/Python semantic-version mismatch.
- `NativeCallError` records `operation`, integer `status`, and the native status
  description.

Platform wheels place the shared library in `metasona/_native`. Developers
can override discovery with `METASONA_LIBRARY` set to an exact compatible
library file. The loader verifies both `ms_abi_version() == 1` and an exact
`ms_version_string()` match with the Python package version.

The loader fails closed when a bundled or explicitly selected library cannot
be loaded. It does not silently substitute a same-named system library. A
controlled source installation may opt in to system lookup with
`METASONA_ALLOW_SYSTEM_LIBRARY=1`; the same ABI and semantic-version checks
still apply. Prefer the exact `METASONA_LIBRARY` path for development.

The environment override is intended for development and controlled
deployments. Treat changes to it like changes to executable search paths; do
not populate it from untrusted input.
