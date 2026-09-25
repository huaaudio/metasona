# SPDX-License-Identifier: Apache-2.0
# MetaSona author: Jiahua Zhang, September 2026.

"""Psychoacoustic metrics with a stable native C core.

Importing this package never loads native code. The bundled shared library is
located and ABI-checked only when a numerical API is first called.
"""

from ._api import (
    loudness_from_levels,
    native_version,
    roughness_daniel_weber,
    sharpness_din45692,
    stationary_loudness,
    third_octave_centres_hz,
    time_varying_loudness,
    tonality_aures,
)
from ._types import (
    LoudnessResult,
    NativeStatus,
    RoughnessResult,
    SharpnessResult,
    SoundField,
    TimeVaryingLoudnessResult,
    TonalityResult,
)
from ._ecma import ecma_tonal_analysis, loudness_ecma, roughness_ecma, tonality_ecma
from ._ecma_types import EcmaLoudnessResult, EcmaRoughnessResult, EcmaTonalAnalysis, EcmaTonalityResult
from .exceptions import (
    MetaSonaError,
    MetaSonaValidationError,
    NativeCallError,
    NativeLibraryError,
)

__all__ = [
    "EcmaLoudnessResult",
    "EcmaRoughnessResult",
    "EcmaTonalAnalysis",
    "EcmaTonalityResult",
    "ecma_tonal_analysis",
    "loudness_ecma",
    "roughness_ecma",
    "tonality_ecma",
    "LoudnessResult",
    "MetaSonaError",
    "MetaSonaValidationError",
    "NativeCallError",
    "NativeLibraryError",
    "NativeStatus",
    "RoughnessResult",
    "SharpnessResult",
    "SoundField",
    "TimeVaryingLoudnessResult",
    "TonalityResult",
    "loudness_from_levels",
    "native_version",
    "roughness_daniel_weber",
    "sharpness_din45692",
    "stationary_loudness",
    "third_octave_centres_hz",
    "time_varying_loudness",
    "tonality_aures",
]

__version__ = "0.2.1"
