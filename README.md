# MetaSona

Psychoacoustic metrics in C, with a typed Python interface.

MetaSona packages compiled C implementations in a Python wheel to reduce
computation time for repeated psychoacoustic analysis, while keeping a simple
NumPy interface. The speed comes from the native implementation; the wheel
makes that implementation easy to install.

| Function | Calculates | Unit |
|---|---|---|
| `stationary_loudness` | Zwicker loudness from audio | sone |
| `loudness_from_levels` | Zwicker loudness from 28 third-octave levels | sone |
| `time_varying_loudness` | Zwicker loudness over time | sone |
| `roughness_daniel_weber` | Daniel–Weber roughness | asper |
| `tonality_aures` | Aures tonality | dimensionless |
| `sharpness_din45692` | Sharpness from specific loudness | acum |
| `loudness_ecma` | ECMA-418-2:2025 loudness, including tonal/noise separation | sone_HMS |
| `roughness_ecma` | ECMA-418-2:2025 roughness | asper_HMS |
| `tonality_ecma` | ECMA-418-2:2025 psychoacoustic tonality | tu_HMS |

The current version is experimental, particularly roughness and tonality. Standards
conformance has not been established.

The metrics draw on [SQAT revision e6228b789fc9](https://github.com/ggrecow/SQAT/tree/e6228b789fc9a22251314f95678b9b1e08e60c55/psychoacoustic_metrics).
Aures tonality uses 250 ms windows and 125 ms hops. MetaSona returns every complete
frame with centre timestamps; see [implementation credits and differences](THIRD_PARTY.md).

## Computation time

Measured September 2026 with the SQAT-aligned MetaSona wheel. Median of five
warmed calls for 10 s of 48 kHz audio (1 kHz carrier, 70 Hz AM plus noise,
60 dB SPL), Ryzen 9 5950X, Windows. MetaSona: MSVC Release; MoSQITo: 1.2.1;
imports and audio I/O are excluded.

| Metric | MetaSona | MoSQITo | Speedup vs MoSQITo |
|---|---:|---:|---:|
| Stationary loudness | 0.090 s | 0.220 s | 2.4× |
| Time-varying loudness | 0.233 s | 17.092 s | 73.5× |
| Roughness | 1.448 s | 22.754 s | 15.7× |
| Aures tonality | 0.572 s | Not available | — |

These are workload timings; models and outputs differ. Native filtering and channel processing contribute to
speed alongside FFTs; compiler and framing differences also affect timings.

## Python

MetaSona requires Python 3.9–3.14. Install the published package from PyPI
with pip:

```sh
python -m pip install metasona
```

Or add it to a uv project:

```sh
uv add metasona
```

Prebuilt wheels are published for Windows x64, Linux x64/ARM64, and macOS
Intel/Apple silicon, so these installations do not require a C compiler.

To install from a source checkout instead, use `python -m pip install .` or
`uv add /path/to/metasona`. Building from source requires a C11 compiler.

Pass **calibrated pressure in pascals**, not uncalibrated audio samples.
The ECMA functions accept mono `(samples,)` or stereo `(samples, 2)` arrays;
the other functions accept mono arrays.
Python accepts integer sample rates from 8 to 192 kHz and resamples to 48 kHz
when needed. Results include their units and read-only NumPy arrays.

Applications that receive calibrated mono audio in chunks can use
`RollingAnalyzer` to compute selected metrics over bounded rolling windows.
The analyzer is synchronous and device-independent, so callers retain control
of capture and concurrency. See the [rolling analysis API](docs/python-api.md#rolling-analysis)
and [complete rolling example](examples/metasona_rolling.py).

For ECMA loudness and tonality together, use
`result = ms.ecma_tonal_analysis(pressure_pa, fs)` and read `result.loudness`
and `result.tonality`. Loudness/tonality require at least 304 ms of audio;
roughness requires 320 ms. Representative values exclude the initial filter
transient. ECMA loudness implements Section 8, whereas MoSQITo's
`loudness_ecma` implements the Section 5 basis calculation. ECMA roughness
does not apply the optional entropy weighting. Function and result-class
docstrings describe output shapes, units, timing and binaural combination.

```python
import numpy as np
import metasona as ms

# One second of a 1 kHz tone at 60 dB SPL (0.02 Pa RMS).
fs = 48_000
t = np.arange(fs) / fs
pressure_pa = np.sqrt(2) * 0.02 * np.sin(2 * np.pi * 1_000 * t)

loudness = ms.stationary_loudness(pressure_pa, fs, sound_field="free")
sharpness = ms.sharpness_din45692(loudness.specific_loudness_sone_per_bark)

print(f"Loudness: {loudness.loudness_sone:.3f} sone")
print(f"Sharpness: {sharpness.sharpness_acum:.3f} acum")
```

See the [Python API](https://github.com/huaaudio/metasona/blob/main/docs/python-api.md)
and [complete example](https://github.com/huaaudio/metasona/blob/main/examples/metasona_example.py).

## C library

Requires CMake 3.21 or newer and a C11 compiler:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
cmake --install build --config Release --prefix install
```

Include `<metasona/metasona.h>` and link the installed CMake target
`MetaSona::metasona`. Native signal functions require 48 kHz input.
Use `-DMS_BUILD_SHARED=OFF` for a static library.

See the [C API](https://github.com/huaaudio/metasona/blob/main/docs/c-api.md)
and [C example](https://github.com/huaaudio/metasona/blob/main/examples/metasona_example.c).

## License and credits

[GPL-3.0](https://github.com/huaaudio/metasona/blob/main/LICENSE) for the combined
distribution, with retained Apache-2.0, BSD-3-Clause and MIT component notices. See
[third-party credits](https://github.com/huaaudio/metasona/blob/main/THIRD_PARTY.md).

## Acknowledgements

Developed by Jiahua Zhang with assistance from OpenAI Codex.

Developed during PhD research within the [METAVISION](https://www.heu-metavision.eu/)
MSCA Doctoral Network.

The European Commission is gratefully acknowledged for their support of the Horizon Europe DN METAVISION project (GA 101072415). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union. The European Union cannot be held responsible for them.
