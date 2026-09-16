# SPDX-License-Identifier: Apache-2.0
# MetaSona author: Jiahua Zhang, September 2026.

"""Validated public Python API backed by the versioned C ABI."""

from __future__ import annotations

import ctypes
from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._abi import DoublePointer, NativeLibrary, get_native_library
from ._resample import NATIVE_SAMPLE_RATE_HZ, to_native_rate
from ._types import (
    LoudnessResult,
    RoughnessResult,
    SharpnessResult,
    SoundField,
    TimeVaryingLoudnessResult,
    TonalityResult,
)
from ._validation import (
    levels_28,
    mono_pressure,
    positive_sample_rate,
    sound_field,
    specific_loudness_240,
    time_skip,
)
from .exceptions import MetaSonaValidationError, NativeCallError

_BARK_BANDS = 240
_FIELD_TO_NATIVE = {SoundField.FREE: 0, SoundField.DIFFUSE: 1}


def _pointer(array: NDArray[np.float64]) -> Any:
    return array.ctypes.data_as(DoublePointer)


def _check(binding: NativeLibrary, operation: str, status: int) -> None:
    if int(status) != 0:
        raise NativeCallError(operation, int(status), binding.status_message(int(status)))


def _prepare_signal(
    pressure_pa: ArrayLike, sample_rate_hz: int
) -> tuple[NDArray[np.float64], float]:
    rate = positive_sample_rate(sample_rate_hz)
    signal = mono_pressure(pressure_pa)
    duration_s = signal.size / rate
    converted = to_native_rate(signal, rate)
    if converted.size == 0 or not np.all(np.isfinite(converted)):
        raise MetaSonaValidationError("resampling did not produce finite pressure samples")
    return converted, duration_s


def _frame_count(
    binding: NativeLibrary,
    operation: str,
    function: Callable[..., int],
    sample_count: int,
) -> int:
    count = ctypes.c_size_t()
    status = function(sample_count, NATIVE_SAMPLE_RATE_HZ, ctypes.byref(count))
    _check(binding, operation, status)
    return int(count.value)


def _field(value: SoundField | str) -> int:
    return _FIELD_TO_NATIVE[sound_field(value)]


def _frame_selection_mask(time_s: NDArray[np.float64], skip_s: float) -> NDArray[np.bool_]:
    mask = time_s >= skip_s
    if not np.any(mask):
        raise MetaSonaValidationError("time_skip_s removes every complete output frame")
    return mask


def native_version() -> str:
    """Return the semantic version reported by the loaded native library."""
    return get_native_library().version


def third_octave_centres_hz() -> NDArray[np.float64]:
    """Return a read-only copy of the 28 nominal third-octave centre frequencies."""
    result = get_native_library().third_octave_centres_hz()
    result.setflags(write=False)
    return result


def loudness_from_levels(
    levels_db_spl: ArrayLike,
    *,
    sound_field: SoundField | str = SoundField.FREE,
) -> LoudnessResult:
    """Compute stationary loudness from 28 third-octave levels.

    ``levels_db_spl`` is ordered from 25 Hz through 12.5 kHz and expressed in
    dB SPL relative to 20 micropascals. Values from 25 through 250 Hz must not
    exceed 120 dB SPL. Input data is copied and never mutated.
    """
    levels = levels_28(levels_db_spl)
    native_field = _field(sound_field)
    binding = get_native_library()
    total = ctypes.c_double()
    specific = np.empty(_BARK_BANDS, dtype=np.float64)
    status = binding.lib.ms_loudness_from_levels(
        _pointer(levels),
        levels.size,
        native_field,
        ctypes.byref(total),
        _pointer(specific),
        specific.size,
    )
    _check(binding, "ms_loudness_from_levels", status)
    return LoudnessResult(float(total.value), specific, binding.bark_axis())


def stationary_loudness(
    pressure_pa: ArrayLike,
    sample_rate_hz: int,
    *,
    sound_field: SoundField | str = SoundField.FREE,
) -> LoudnessResult:
    """Compute ISO-532-1:2017-targeted stationary loudness.

    Audio must be calibrated mono acoustic pressure in pascals. The Python
    boundary accepts practical integer audio rates from 8 to 192 kHz and
    converts them to exactly 48 kHz using the documented polyphase FIR
    boundary; the native loudness kernel only sees 48 kHz.
    """
    signal, _ = _prepare_signal(pressure_pa, sample_rate_hz)
    native_field = _field(sound_field)
    binding = get_native_library()
    total = ctypes.c_double()
    specific = np.empty(_BARK_BANDS, dtype=np.float64)
    status = binding.lib.ms_loudness_stationary(
        _pointer(signal),
        signal.size,
        NATIVE_SAMPLE_RATE_HZ,
        native_field,
        ctypes.byref(total),
        _pointer(specific),
        specific.size,
    )
    _check(binding, "ms_loudness_stationary", status)
    return LoudnessResult(float(total.value), specific, binding.bark_axis())


def time_varying_loudness(
    pressure_pa: ArrayLike,
    sample_rate_hz: int,
    *,
    sound_field: SoundField | str = SoundField.FREE,
    time_skip_s: float = 0.0,
) -> TimeVaryingLoudnessResult:
    """Compute 2 ms time-varying loudness frames.

    Timestamps are frame centres: the first 2 ms frame is labelled 1 ms.
    ``time_skip_s`` filters complete results and must leave at least one frame.
    """
    signal, duration_s = _prepare_signal(pressure_pa, sample_rate_hz)
    skip_s = time_skip(time_skip_s, duration_s)
    native_field = _field(sound_field)
    binding = get_native_library()
    frame_count = _frame_count(
        binding,
        "ms_loudness_time_frame_count",
        binding.lib.ms_loudness_time_frame_count,
        signal.size,
    )
    total = np.empty(frame_count, dtype=np.float64)
    specific = np.empty((frame_count, _BARK_BANDS), dtype=np.float64)
    written = ctypes.c_size_t()
    status = binding.lib.ms_loudness_time(
        _pointer(signal),
        signal.size,
        NATIVE_SAMPLE_RATE_HZ,
        native_field,
        _pointer(total),
        total.size,
        _pointer(specific.reshape(-1)),
        specific.size,
        ctypes.byref(written),
    )
    _check(binding, "ms_loudness_time", status)
    if written.value != frame_count:
        raise NativeCallError(
            "ms_loudness_time",
            -1,
            f"frame-count contract violation: queried {frame_count}, wrote {written.value}",
        )
    centres = (np.arange(frame_count, dtype=np.float64) + 0.5) * 0.002
    selected = _frame_selection_mask(centres, skip_s)
    return TimeVaryingLoudnessResult(
        centres[selected],
        total[selected],
        specific[selected, :],
        binding.bark_axis(),
    )


def roughness_daniel_weber(
    pressure_pa: ArrayLike,
    sample_rate_hz: int,
    *,
    time_skip_s: float = 0.0,
) -> RoughnessResult:
    """Compute the experimental Daniel--Weber-targeted tabulated model in asper.

    Frames are 200 ms with 50 percent overlap. Timestamps are frame centres,
    starting at 100 ms. Uses the SQAT-aligned psychohelperc kernel and
    all complete frames. Standards validation is pending.
    """
    signal, duration_s = _prepare_signal(pressure_pa, sample_rate_hz)
    skip_s = time_skip(time_skip_s, duration_s)
    binding = get_native_library()
    frame_count = _frame_count(
        binding,
        "ms_roughness_frame_count",
        binding.lib.ms_roughness_frame_count,
        signal.size,
    )
    values = np.empty(frame_count, dtype=np.float64)
    written = ctypes.c_size_t()
    status = binding.lib.ms_roughness_dw(
        _pointer(signal),
        signal.size,
        NATIVE_SAMPLE_RATE_HZ,
        _pointer(values),
        values.size,
        ctypes.byref(written),
    )
    _check(binding, "ms_roughness_dw", status)
    if written.value != frame_count:
        raise NativeCallError(
            "ms_roughness_dw",
            -1,
            f"frame-count contract violation: queried {frame_count}, wrote {written.value}",
        )
    centres = 0.1 + np.arange(frame_count, dtype=np.float64) * 0.1
    selected = _frame_selection_mask(centres, skip_s)
    return RoughnessResult(centres[selected], values[selected])


def tonality_aures(
    pressure_pa: ArrayLike,
    sample_rate_hz: int,
    *,
    sound_field: SoundField | str = SoundField.FREE,
    time_skip_s: float = 0.0,
) -> TonalityResult:
    """Compute experimental Aures-1985-targeted tonality.

    Frames are 250 ms with 50 percent overlap. Timestamps are frame centres,
    starting at 125 ms. Uses the SQAT-aligned psychohelperc kernel and
    all complete frames.
    """
    signal, duration_s = _prepare_signal(pressure_pa, sample_rate_hz)
    skip_s = time_skip(time_skip_s, duration_s)
    native_field = _field(sound_field)
    binding = get_native_library()
    frame_count = _frame_count(
        binding,
        "ms_tonality_frame_count",
        binding.lib.ms_tonality_frame_count,
        signal.size,
    )
    values = np.empty(frame_count, dtype=np.float64)
    written = ctypes.c_size_t()
    status = binding.lib.ms_tonality_aures(
        _pointer(signal),
        signal.size,
        NATIVE_SAMPLE_RATE_HZ,
        native_field,
        _pointer(values),
        values.size,
        ctypes.byref(written),
    )
    _check(binding, "ms_tonality_aures", status)
    if written.value != frame_count:
        raise NativeCallError(
            "ms_tonality_aures",
            -1,
            f"frame-count contract violation: queried {frame_count}, wrote {written.value}",
        )
    centres = 0.125 + np.arange(frame_count, dtype=np.float64) * 0.125
    selected = _frame_selection_mask(centres, skip_s)
    return TonalityResult(centres[selected], values[selected])


def sharpness_din45692(specific_loudness_sone_per_bark: ArrayLike) -> SharpnessResult:
    """Compute DIN 45692 sharpness from one 240-point loudness pattern."""
    specific = specific_loudness_240(specific_loudness_sone_per_bark)
    binding = get_native_library()
    result = ctypes.c_double()
    status = binding.lib.ms_sharpness_din(_pointer(specific), specific.size, ctypes.byref(result))
    _check(binding, "ms_sharpness_din", status)
    return SharpnessResult(float(result.value))
