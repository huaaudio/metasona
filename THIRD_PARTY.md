# Third-party credits

The combined distribution is licensed under [GPL-3.0](LICENSE), following the
current SQAT tonality adaptation. MetaSona's original code retains its
[Apache-2.0](LICENSES/Apache-2.0.txt) notice. BSD-3-Clause and MIT components
retain their notices in [NOTICE](NOTICE) and `LICENSES/`.

## PocketFFT

The C backend in `src/third_party/pocketfft/` is vendored unchanged from
[mreineck/pocketfft](https://github.com/mreineck/pocketfft/tree/81d171a6d5562e3aaa2c73489b70f564c633ff81),
commit `81d171a6d5562e3aaa2c73489b70f564c633ff81` (C/master branch).
Copyright (C) 2010-2019 Max-Planck-Society; author Martin Reinecke.
It is distributed under the [BSD-3-Clause license](LICENSES/BSD-3-Clause-PocketFFT.txt).
The upstream notice is included both beside the sources and in the
distribution's license directory. These C files match the supplied
psychohelperc copies after normalizing line endings.

The reused kernels call PocketFFT directly and own their plans per call.
Build definitions prefix backend symbols and hide them from the shared ABI;
the upstream source is not edited. No global plan cache or OpenMP runtime is
introduced. See the pinned hashes in `src/third_party/pocketfft/README.md`.

## MoSQITo

[MoSQITo](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/README.md), revision `d990c33f94f1b2050b2d811dea492219ebe97b30`,
is distributed under [Apache-2.0](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/LICENSE).
Credit: Eomys and the MoSQITo contributors.

Model data correspond to these upstream files:

- [Loudness filter coefficients](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwtv/_third_octave_levels.py): `src/ms_loudness_data.h`.
- [Core loudness](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwst/_main_loudness.py) and [upper slopes](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwst/_calc_slopes.py): tables in `src/ms_loudness.c`.
- [Bark conversion](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/utils/conversion/freq2bark.py), [hearing threshold](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/utils/LTQ.py), [ear transfer](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/roughness/roughness_dw/_ear_filter_coeff.py) and [channel weighting](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/roughness/roughness_dw/_gzi_weighting.py): corresponding data in the reused `src/roughness_dw.c`.

Tables were checked element by element. One deliberate loudness difference is
retained: `MS_UPPER_SLOPES[12][4]` is `0.22`, whereas this MoSQITo revision
uses `0.24`. This value comes from the supplied psychohelperc/ISO_532-1.c table; it is not a fixture adjustment.
MetaSona reorganizes data into const C arrays and uses its own C API,
allocation, FFT and error-handling infrastructure. MoSQITo is not a runtime
dependency.

## SQAT

The ECMA-418-2:2025 kernels (`src/ms_ecma*.c`) adapt the SQAT/RefMap
implementations of loudness, roughness, tonality and shared hearing-model
utilities by Mike JB Lotinga and Matt Torjussen, University of Salford,
RefMap project. These adaptations retain GPL-3.0. They use call-local
buffers, portable C FFTs and checked C/Python interfaces. MoSQITo's
Apache-2.0 Section 5 implementation also informed the shared primitives.
Detailed validation and reference comparisons are maintained outside the
package repository in the development workspace.

The roughness port uses the standard's exact 20 ms grid and nominal
calibration factor 0.0180685; SQAT uses a grid ending at the signal duration
and a slightly adjusted calibration. Equal-prominence peaks are selected
deterministically by lower frequency. The optional entropy weighting is
not applied. Tonal frequency is selected from the positive-frequency half
of the real autocorrelation spectrum; this resolves conjugate peak ties.
The loudness combination uses a scaled power sum to avoid underflow near
silence. The autocorrelation explicitly zeros overlaps with zero energy
and enforces the mathematical [0,1] bound after rectification, preventing
FFT roundoff from creating false components during startup. These are
numerical/implementation differences, not new models.

The following files in [SQAT](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/README.md), revision
`e6228b789fc9a22251314f95678b9b1e08e60c55`, carry explicit file-level licenses:

| Upstream implementation | License | Use in MetaSona |
|---|---|---|
| [Loudness](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Loudness_ISO532_1/Loudness_ISO532_1.m) and [third-octave filtering](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/sound_level_meter/Do_OB13_ISO532_1.m) | [BSD-3-Clause](LICENSES/BSD-3-Clause-SQAT.txt), Ella Manor | Licensed references for filter-bank and temporal loudness processing |
| [Daniel–Weber roughness](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Roughness_Daniel1997/Roughness_Daniel1997.m) and [modulation weights](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Roughness_Daniel1997/private/Get_Hweight_roughness.m) | [MIT](LICENSES/MIT-SQAT.txt), Dik Hermes | Roughness processing reference and H2/H5/H16/H21/H42 table data |

The roughness port also follows the explicitly MIT-licensed SQAT utilities
`Terhardt_filterbank.m`, `Terhardt_filterbank_params.m`, `Get_Bark.m`, and
`calculate_a0.m`. Contributors include Alejandro Osses, Sergio Aguirre and
Gil Felix Greco. MetaSona does not bundle GPL-licensed SQAT wrappers.
The licensed upstream implementations above provide explicit attribution for
the corresponding processing and data. The user-authored
`psychohelperc/roughness_dw.c` and `.h` are reused directly as
`src/roughness_dw.c` and `.h`, with the author's permission. MetaSona adds
checked cleanup, FFT/numerical error propagation and the final complete frame.
`src/ms_signal_metrics.c` provides validation and atomic output updates.
This is source-informed code, not a claim of independent or clean-room origin.

Roughness follows this pinned SQAT revision's periodic Blackman window,
frequency-bin mapping, Terhardt excitation, modulation filters and calibration
`0.25 / 1.006197`. The latter is SQAT's measured reference normalization,
not an extra MetaSona correction. The public stationary-loudness function
continues to use the whole record.

## PA Aures tonality

The author's `psychohelperc/tonality_aures1985.c` and `.h` are reused directly
as `src/tonality_aures1985.c` and `.h`, with the author's permission. They are
derived from the supplied `PA/Tonality_Aures1985/Tonality_Aures1985.m`, by
Gil Felix Greco (2020, updated 2023 and 2025). Its file-level BSD-3-Clause
notice is retained in `LICENSES/BSD-3-Clause-PA-Tonality.txt`.

MetaSona adapts the loudness callback to its existing core, includes the final
complete frame and checks FFT/numerical failures. The supplied C port's
phase-preserving replacement and per-tone masking are retained. The port now
follows [SQAT's pinned tonality implementation](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Tonality_Aures1985/Tonality_Aures1985.m):
250 ms Hann windows, sinusoidal and narrowband-region detection, residual-noise
masking, window-width correction, windowed loudness and calibration `1.1055`.
The narrowband implementation is in `src/tonality_narrowband.h`. Credit also
goes to the September 2026 SQAT contributors, Sergio Aguirre and Gil Felix Greco.
Unlike the older PA source, this revision explicitly licenses its tonality
file under GPL-3.0. The translated detector and updated kernel therefore carry
GPL-3.0 notices; the combined MetaSona distribution uses GPL-3.0. SQAT discloses
AI assistance by Claude Fable 5.1 and Opus 5 for its September 2026 changes.

### Remaining implementation differences

- MetaSona returns all complete frames with centre timestamps; SQAT's wrappers
  have different trailing-frame and timestamp conventions.
- Python resamples to 48 kHz. Native kernels require 48 kHz.
- Tonality uses MetaSona's existing ISO loudness core, skipping the first 5%
  of each window when averaging filter energy. Exact MATLAB agreement is not
  established. Source-based Python comparisons are development checks, not an
  executed MATLAB/Octave validation.
- Invalid bandwidths, overflow and FFT failures produce errors rather than
  SQAT's substituted 1 Hz bandwidth or a plausible zero result.
- Roughness retains SQAT's zero clamp on non-decaying upper excitation slopes.
  It can occur above approximately 121 dB per component at 1 kHz. The C API
  has no warning channel for this extrapolation beyond the model's usual range.
- FFT rounding can affect correlation of almost constant channel envelopes;
  numerical agreement does not imply bitwise equality or standards conformance.
