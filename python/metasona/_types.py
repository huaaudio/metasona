# SPDX-License-Identifier: Apache-2.0
# MetaSona authors: Jiahua Zhang and Codex, September 2026.

"""Public enums and immutable result containers."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, IntEnum
from typing import TypeAlias

import numpy as np
from numpy.typing import NDArray

FloatArray: TypeAlias = NDArray[np.float64]


class SoundField(str, Enum):
    """Sound-field correction used by loudness-dependent calculations."""

    FREE = "free"
    DIFFUSE = "diffuse"


class NativeStatus(IntEnum):
    """Stable status values from the C ABI."""

    OK = 0
    NULL_POINTER = 1
    INVALID_ARGUMENT = 2
    INVALID_SAMPLE_RATE = 3
    INVALID_SOUND_FIELD = 4
    NONFINITE_INPUT = 5
    INPUT_TOO_SHORT = 6
    OUTPUT_TOO_SMALL = 7
    SIZE_OVERFLOW = 8
    ALLOCATION = 9
    NUMERICAL = 10


def _immutable_float_array(value: FloatArray) -> FloatArray:
    result = np.array(value, dtype=np.float64, copy=True, order="C")
    result.setflags(write=False)
    return result


@dataclass(frozen=True, slots=True)
class LoudnessResult:
    """Stationary loudness and its 0.1-to-24 Bark specific pattern."""

    loudness_sone: float
    specific_loudness_sone_per_bark: FloatArray
    bark_axis: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "specific_loudness_sone_per_bark",
            _immutable_float_array(self.specific_loudness_sone_per_bark),
        )
        object.__setattr__(self, "bark_axis", _immutable_float_array(self.bark_axis))


@dataclass(frozen=True, slots=True)
class TimeVaryingLoudnessResult:
    """Frame-centred time-varying total and specific loudness."""

    time_s: FloatArray
    loudness_sone: FloatArray
    specific_loudness_sone_per_bark: FloatArray
    bark_axis: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_s", _immutable_float_array(self.time_s))
        object.__setattr__(self, "loudness_sone", _immutable_float_array(self.loudness_sone))
        object.__setattr__(
            self,
            "specific_loudness_sone_per_bark",
            _immutable_float_array(self.specific_loudness_sone_per_bark),
        )
        object.__setattr__(self, "bark_axis", _immutable_float_array(self.bark_axis))


@dataclass(frozen=True, slots=True)
class RoughnessResult:
    """Frame-centred Daniel--Weber roughness time series."""

    time_s: FloatArray
    roughness_asper: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_s", _immutable_float_array(self.time_s))
        object.__setattr__(self, "roughness_asper", _immutable_float_array(self.roughness_asper))


@dataclass(frozen=True, slots=True)
class TonalityResult:
    """Frame-centred Aures-1985-targeted tonality time series."""

    time_s: FloatArray
    tonality: FloatArray

    def __post_init__(self) -> None:
        object.__setattr__(self, "time_s", _immutable_float_array(self.time_s))
        object.__setattr__(self, "tonality", _immutable_float_array(self.tonality))


@dataclass(frozen=True, slots=True)
class SharpnessResult:
    """DIN 45692 sharpness in acum."""

    sharpness_acum: float
