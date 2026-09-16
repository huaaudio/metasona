# SPDX-License-Identifier: Apache-2.0
# MetaSona authors: Jiahua Zhang and Codex, September 2026.

"""Validation shared by all Python API entry points."""

from __future__ import annotations

import math
from numbers import Integral, Real
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._types import SoundField
from .exceptions import MetaSonaValidationError

UINT32_MAX = 2**32 - 1
MIN_SAMPLE_RATE_HZ = 8_000
MAX_SAMPLE_RATE_HZ = 192_000


def positive_sample_rate(value: Any) -> int:
    """Return a practical audio sample rate suitable for the uint32 C ABI."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral):
        raise MetaSonaValidationError("sample_rate_hz must be a positive integer")
    result = int(value)
    if result < MIN_SAMPLE_RATE_HZ or result > MAX_SAMPLE_RATE_HZ:
        raise MetaSonaValidationError(
            f"sample_rate_hz must be between {MIN_SAMPLE_RATE_HZ} and {MAX_SAMPLE_RATE_HZ}"
        )
    return result


def _real_numeric_array(value: ArrayLike, name: str) -> NDArray[np.float64]:
    """Convert ordinary real numeric arrays without silently coercing categories."""
    try:
        source = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise MetaSonaValidationError(f"{name} must be a real numeric array") from exc
    if source.dtype.kind in {"b", "c", "O", "S", "U", "V"}:
        raise MetaSonaValidationError(
            f"{name} must contain real numeric values, not boolean, complex, object, or text data"
        )
    try:
        return np.asarray(source, dtype=np.float64)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MetaSonaValidationError(f"{name} must be a real numeric array") from exc


def mono_pressure(audio: ArrayLike) -> NDArray[np.float64]:
    """Validate finite, nonempty mono pressure samples and copy to C order."""
    array = _real_numeric_array(audio, "pressure_pa")
    if array.ndim != 1:
        raise MetaSonaValidationError(
            "pressure_pa must be one-dimensional mono audio; select or mix channels explicitly"
        )
    if array.size == 0:
        raise MetaSonaValidationError("pressure_pa must not be empty")
    if not np.all(np.isfinite(array)):
        raise MetaSonaValidationError("pressure_pa must contain only finite values")
    return np.array(array, dtype=np.float64, copy=True, order="C")


def levels_28(levels_db_spl: ArrayLike) -> NDArray[np.float64]:
    """Validate the native 28-band third-octave level representation."""
    array = _real_numeric_array(levels_db_spl, "levels_db_spl")
    if array.shape != (28,):
        raise MetaSonaValidationError(
            "levels_db_spl must contain exactly 28 one-third-octave levels"
        )
    if not np.all(np.isfinite(array)):
        raise MetaSonaValidationError("levels_db_spl must contain only finite values")
    if np.any(array[:11] > 120.0):
        raise MetaSonaValidationError(
            "levels_db_spl values from 25 through 250 Hz must not exceed 120 dB SPL"
        )
    return np.array(array, dtype=np.float64, copy=True, order="C")


def specific_loudness_240(values: ArrayLike) -> NDArray[np.float64]:
    """Validate one nonnegative 240-point specific-loudness pattern."""
    array = _real_numeric_array(values, "specific_loudness")
    if array.shape != (240,):
        raise MetaSonaValidationError(
            "specific_loudness must contain exactly 240 values from 0.1 to 24 Bark"
        )
    if not np.all(np.isfinite(array)):
        raise MetaSonaValidationError("specific_loudness must contain only finite values")
    if np.any(array < 0.0):
        raise MetaSonaValidationError("specific_loudness must not contain negative values")
    return np.array(array, dtype=np.float64, copy=True, order="C")


def sound_field(value: SoundField | str) -> SoundField:
    """Parse a strict public sound-field value."""
    if isinstance(value, SoundField):
        return value
    if isinstance(value, str):
        try:
            return SoundField(value)
        except ValueError as exc:
            raise MetaSonaValidationError("sound_field must be 'free' or 'diffuse'") from exc
    raise MetaSonaValidationError("sound_field must be 'free' or 'diffuse'")


def time_skip(value: Any, duration_s: float) -> float:
    """Validate a nonnegative finite skip that leaves some input audio."""
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise MetaSonaValidationError("time_skip_s must be a finite nonnegative number")
    result = float(value)
    if not math.isfinite(result) or result < 0.0:
        raise MetaSonaValidationError("time_skip_s must be a finite nonnegative number")
    if result >= duration_s:
        raise MetaSonaValidationError("time_skip_s must be shorter than the input duration")
    return result
