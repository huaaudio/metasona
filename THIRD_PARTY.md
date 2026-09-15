# Third-party credits

MetaSona's original code is licensed under [Apache-2.0](LICENSE).
The distribution also contains adaptations covered by BSD-3-Clause and MIT;
their notices are included in [NOTICE](NOTICE) and `LICENSES/`.

## MoSQITo

[MoSQITo](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/README.md), revision `d990c33f94f1b2050b2d811dea492219ebe97b30`,
is distributed under [Apache-2.0](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/LICENSE).
Credit: Eomys and the MoSQITo contributors.

Model data correspond to these upstream files:

- [Loudness filter coefficients](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwtv/_third_octave_levels.py): `src/ms_loudness_data.h`.
- [Core loudness](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwst/_main_loudness.py) and [upper slopes](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/loudness/loudness_zwst/_calc_slopes.py): tables in `src/ms_loudness.c`.
- [Bark conversion](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/utils/conversion/freq2bark.py), [hearing threshold](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/utils/LTQ.py), [ear transfer](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/roughness/roughness_dw/_ear_filter_coeff.py) and [channel weighting](https://github.com/Eomys/MoSQITo/blob/d990c33f94f1b2050b2d811dea492219ebe97b30/mosqito/sq_metrics/roughness/roughness_dw/_gzi_weighting.py): tables in `src/ms_roughness_data.h`.

Tables were checked element by element. One deliberate loudness difference is
retained: `MS_UPPER_SLOPES[12][4]` is `0.22`, whereas this MoSQITo revision
uses `0.24`. This preserves the earlier reference comparison behavior.
MetaSona reorganizes data into const C arrays and uses its own C API,
allocation, FFT and error-handling infrastructure. MoSQITo is not a runtime
dependency.

## SQAT

The following files in [SQAT](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/README.md), revision
`e6228b789fc9a22251314f95678b9b1e08e60c55`, carry explicit file-level licenses:

| Upstream implementation | License | Use in MetaSona |
|---|---|---|
| [Loudness](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Loudness_ISO532_1/Loudness_ISO532_1.m) and [third-octave filtering](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/sound_level_meter/Do_OB13_ISO532_1.m) | [BSD-3-Clause](LICENSES/BSD-3-Clause-SQAT.txt), Ella Manor | Licensed references for filter-bank and temporal loudness processing |
| [Daniel–Weber roughness](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Roughness_Daniel1997/Roughness_Daniel1997.m) and [modulation weights](https://github.com/ggrecow/SQAT/blob/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics/Roughness_Daniel1997/private/Get_Hweight_roughness.m) | [MIT](LICENSES/MIT-SQAT.txt), Dik Hermes | Roughness processing reference and H2/H5/H16/H21/H42 table data |

MetaSona does not include SQAT's GPL-licensed file wrappers or utilities.
The supplied legacy C informed earlier development; the licensed upstream
implementations above now provide explicit attribution for the corresponding
processing and data. The original C program is not bundled. This is a
source-informed implementation, not a claim of independent or clean-room origin.

Behavior differs from current SQAT: stationary analysis uses the whole record;
roughness uses a symmetric Blackman window, includes the last complete frame,
retains the 500 Hz grid cutoff for H5/H16/H21/H42, and uses factor `0.25` with
MoSQITo channel weights. These choices do not establish numerical equivalence
or standards conformance. See [model conventions](docs/algorithm-notes.md).

## Python dependencies

NumPy 2.2.6 and SciPy 1.15.3 are installed separately and retain their own BSD
licenses and bundled-component notices. They provide arrays and resampling;
their source and binaries are not vendored in this repository. Build tools
are also supplied separately.
