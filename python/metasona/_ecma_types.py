# SPDX-License-Identifier: Apache-2.0
# MetaSona author: Jiahua Zhang, 2026.
"""Typed ECMA-418-2:2025 result containers."""
from __future__ import annotations

from dataclasses import dataclass, fields
from typing import Any

import numpy as np

from ._types import FloatArray, _immutable_float_array


def _freeze_arrays(result: Any) -> None:
    for field in fields(result):
        value = getattr(result, field.name)
        if isinstance(value, np.ndarray):
            object.__setattr__(result, field.name, _immutable_float_array(value))


@dataclass(frozen=True, slots=True)
class EcmaLoudnessResult:
    """Section 8 loudness in sone_HMS, with the 53-band Bark_HMS pattern.

    Mono arrays have shape (frames,) or (frames, 53). Stereo adds a final
    axis of length two (left, right). Representative values are scalars for
    mono and (2,) arrays for stereo. Binaural arrays are present only for
    stereo and combine each band's left/right values by quadratic mean.
    ``mean_loudness_sone`` is the power mean after the first 304 ms.
    """

    time_s: FloatArray
    loudness_sone: FloatArray
    specific_loudness_sone_per_bark: FloatArray
    mean_loudness_sone: float | FloatArray
    mean_specific_loudness_sone_per_bark: FloatArray
    tonal_loudness_sone_per_bark: FloatArray
    noise_loudness_sone_per_bark: FloatArray
    bark_axis: FloatArray
    centre_frequencies_hz: FloatArray
    binaural_loudness_sone: FloatArray | None = None
    binaural_specific_loudness_sone_per_bark: FloatArray | None = None
    binaural_mean_loudness_sone: float | None = None

    def __post_init__(self) -> None:
        _freeze_arrays(self)


@dataclass(frozen=True, slots=True)
class EcmaTonalityResult:
    """Section 6 tonality in tu_HMS, evaluated separately for each ear.

    Shape conventions match EcmaLoudnessResult. Specific tonality has units
    tu_HMS (not tu_HMS/Bark); total tonality is the band maximum. The mean
    excludes the first 304 ms and values <= 0.02 tu_HMS. Mean frequencies
    use the corresponding active frames. Frequency is zero when no tonal
    component is present. No binaural tonality combination is defined here.
    """

    time_s: FloatArray
    tonality_tu: FloatArray
    specific_tonality_tu: FloatArray
    mean_tonality_tu: float | FloatArray
    mean_specific_tonality_tu: FloatArray
    tonal_frequency_hz: FloatArray
    specific_tonal_frequency_hz: FloatArray
    mean_specific_tonal_frequency_hz: FloatArray
    tonal_loudness_sone_per_bark: FloatArray
    noise_loudness_sone_per_bark: FloatArray
    bark_axis: FloatArray
    centre_frequencies_hz: FloatArray

    def __post_init__(self) -> None:
        _freeze_arrays(self)


@dataclass(frozen=True, slots=True)
class EcmaRoughnessResult:
    """Section 7 roughness in asper_HMS, without optional entropy weighting.

    Shape conventions match EcmaLoudnessResult. Output is sampled at 50 Hz.
    The 90th percentile and average specific pattern exclude the first
    320 ms. Percentiles use linear interpolation between sample probability
    positions (i + 0.5) / n (Hazen convention, as in MATLAB prctile).
    """

    time_s: FloatArray
    roughness_asper: FloatArray
    specific_roughness_asper_per_bark: FloatArray
    roughness_90_asper: float | FloatArray
    mean_specific_roughness_asper_per_bark: FloatArray
    bark_axis: FloatArray
    centre_frequencies_hz: FloatArray
    binaural_roughness_asper: FloatArray | None = None
    binaural_specific_roughness_asper_per_bark: FloatArray | None = None
    binaural_roughness_90_asper: float | None = None

    def __post_init__(self) -> None:
        _freeze_arrays(self)


@dataclass(frozen=True, slots=True)
class EcmaTonalAnalysis:
    """Loudness and tonality calculated together with one native analysis."""

    loudness: EcmaLoudnessResult
    tonality: EcmaTonalityResult
