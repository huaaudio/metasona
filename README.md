# MetaSona

Psychoacoustic metrics in C11, with a typed Python interface.

| Function | Calculates | Unit |
|---|---|---|
| `stationary_loudness` | Zwicker loudness from audio | sone |
| `loudness_from_levels` | Zwicker loudness from 28 third-octave levels | sone |
| `time_varying_loudness` | Zwicker loudness over time | sone |
| `roughness_daniel_weber` | Daniel–Weber roughness | asper |
| `tonality_aures` | Aures tonality | dimensionless |
| `sharpness_din45692` | Sharpness from specific loudness | acum |

Version 0.1.0 is experimental, particularly roughness and tonality. Standards
conformance has not been established.

## Python

Requires Python 3.10–3.13 and a C11 compiler. From this checkout:

```sh
python -m pip install .
```

Or add the checkout to a uv project with `uv add /path/to/metasona`.

Pass **calibrated mono pressure in pascals**, not uncalibrated audio samples.
Python accepts integer sample rates from 8 to 192 kHz and resamples to 48 kHz
when needed. Results include their units and read-only NumPy arrays.

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

See the [C API](https://github.com/huaaudio/metasona/blob/main/docs/c-api.md),
[C example](https://github.com/huaaudio/metasona/blob/main/examples/metasona_example.c),
and [model conventions](https://github.com/huaaudio/metasona/blob/main/docs/algorithm-notes.md).

## License and credits

[Apache-2.0](https://github.com/huaaudio/metasona/blob/main/LICENSE), with BSD-3-Clause
and MIT components from MoSQITo and SQAT. See
[third-party credits](https://github.com/huaaudio/metasona/blob/main/THIRD_PARTY.md).

## Acknowledgements

Developed during PhD research within the [METAVISION](https://www.heu-metavision.eu/)
MSCA Doctoral Network.

The European Commission is gratefully acknowledged for their support of the Horizon Europe DN METAVISION project (GA 101072415). Views and opinions expressed are however those of the authors only and do not necessarily reflect those of the European Union. The European Union cannot be held responsible for them.
