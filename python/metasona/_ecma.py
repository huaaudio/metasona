# SPDX-License-Identifier: GPL-3.0-only
# MetaSona author: Jiahua Zhang, 2026. See NOTICE for reference attribution.
"""ECMA-418-2:2025 interfaces; all signal-analysis kernels execute in C."""
from __future__ import annotations

import ctypes

import numpy as np
from numpy.typing import ArrayLike

from ._abi import get_native_library
from ._api import _check, _field, _frame_count, _pointer
from ._ecma_types import (
    EcmaLoudnessResult, EcmaRoughnessResult, EcmaTonalAnalysis, EcmaTonalityResult,
)
from ._resample import NATIVE_SAMPLE_RATE_HZ, to_native_rate
from ._types import FloatArray, SoundField
from ._validation import _real_numeric_array, positive_sample_rate
from .exceptions import MetaSonaValidationError

_BARK = np.arange(.5, 27, .5)
_CENTRES = (81.9289 / .1618) * np.sinh(.1618 * _BARK)


def _signals(pressure_pa: ArrayLike, sample_rate_hz: int, minimum: int) -> tuple[list[FloatArray], bool]:
    rate = positive_sample_rate(sample_rate_hz)
    source = _real_numeric_array(pressure_pa, "pressure_pa")
    mono = source.ndim == 1
    if not mono and (source.ndim != 2 or source.shape[1] != 2):
        raise MetaSonaValidationError("pressure_pa must have shape (samples,) or (samples, 2)")
    if source.shape[0] * NATIVE_SAMPLE_RATE_HZ < minimum * rate:
        raise MetaSonaValidationError(f"pressure_pa must contain at least {minimum / 48:g} ms")
    if not np.all(np.isfinite(source)):
        raise MetaSonaValidationError("pressure_pa must contain only finite values")
    channels = [source] if mono else [source[:, 0], source[:, 1]]
    result = []
    for channel in channels:
        converted = to_native_rate(np.array(channel, dtype=np.float64, copy=True, order="C"), rate)
        if not np.all(np.isfinite(converted)):
            raise MetaSonaValidationError("resampling did not produce finite pressure samples")
        result.append(np.ascontiguousarray(converted))
    return result, mono


def _scalar_or_array(value: FloatArray) -> float | FloatArray:
    return float(value) if np.ndim(value) == 0 else value


def _power_mean(values: FloatArray) -> FloatArray:
    # Scaling preserves accuracy at very low and high finite values.
    scale = np.max(values, axis=0)
    norm = np.divide(values, scale, out=np.zeros_like(values), where=scale > 0)
    return scale * np.mean(norm ** (1 / np.log10(2)), axis=0) ** np.log10(2)


def _active_mean(values: FloatArray, active: np.ndarray) -> FloatArray:
    count = np.sum(active, axis=0)
    total = np.sum(np.where(active, values, 0), axis=0)
    return np.divide(total, count, out=np.zeros_like(total), where=count > 0)


def ecma_tonal_analysis(
    pressure_pa: ArrayLike, sample_rate_hz: int, *,
    sound_field: SoundField | str = SoundField.FREE,
) -> EcmaTonalAnalysis:
    """Calculate ECMA-418-2:2025 loudness and tonality together.

    Supply calibrated pressure in pascals as mono (samples,) or stereo
    (samples, 2), left then right. Input must be at least 304 ms long.
    Integer audio rates from 8 to 192 kHz are converted to 48 kHz using
    the same polyphase FIR boundary as the other MetaSona metrics.

    This is full Section 8 loudness, including tonal/noise separation;
    MoSQITo's loudness_ecma calculates Section 5 basis loudness instead.
    Output times start at zero, separated by 256/48000 s. The final frame
    may extend less than one output interval past the input duration.
    Representative values discard frames 0..56 (the first 304 ms).
    Calling this function avoids duplicate analysis when both metrics
    are needed. Input is never modified; returned arrays are read-only.
    """
    field = _field(sound_field)
    channels, mono = _signals(pressure_pa, sample_rate_hz, 14592)
    binding = get_native_library()
    frames = _frame_count(binding, "ms_ecma_tonal_frame_count",
                          binding.lib.ms_ecma_tonal_frame_count, channels[0].size)
    planes = []
    for channel in channels:
        result = np.empty((5, frames, 53), dtype=np.float64)
        written = ctypes.c_size_t()
        status = binding.lib.ms_ecma_tonal_analysis(
            _pointer(channel), channel.size, NATIVE_SAMPLE_RATE_HZ, field,
            _pointer(result), result.size, ctypes.byref(written))
        _check(binding, "ms_ecma_tonal_analysis", status)
        planes.append(result)
    data = planes[0] if mono else np.stack(planes, axis=-1)
    tone, noise, freq, spec_tonality, specific = data
    time = np.arange(frames, dtype=np.float64) / 187.5
    total = .5 * specific.sum(axis=1)
    binaural = None if mono else np.hypot(specific[..., 0], specific[..., 1]) / np.sqrt(2)
    bin_total = None if binaural is None else .5 * binaural.sum(axis=1)
    loudness = EcmaLoudnessResult(
        time_s=time, loudness_sone=total, specific_loudness_sone_per_bark=specific,
        mean_loudness_sone=_scalar_or_array(_power_mean(total[57:])),
        mean_specific_loudness_sone_per_bark=_power_mean(specific[57:]),
        tonal_loudness_sone_per_bark=tone, noise_loudness_sone_per_bark=noise,
        bark_axis=_BARK, centre_frequencies_hz=_CENTRES,
        binaural_loudness_sone=bin_total, binaural_specific_loudness_sone_per_bark=binaural,
        binaural_mean_loudness_sone=None if bin_total is None else float(_power_mean(bin_total[57:])),
    )
    total_tonality = spec_tonality.max(axis=1)
    # Frequencies of zero-energy tonal components are not physically defined.
    freq = np.where(tone > 0, freq, 0)
    strongest = spec_tonality.argmax(axis=1)
    total_freq = np.take_along_axis(freq, np.expand_dims(strongest, axis=1), axis=1).squeeze(axis=1)
    total_freq = np.where(total_tonality > 0, total_freq, 0)
    active = spec_tonality[57:] > .02
    tonality = EcmaTonalityResult(
        time_s=time, tonality_tu=total_tonality, specific_tonality_tu=spec_tonality,
        mean_tonality_tu=_scalar_or_array(_active_mean(total_tonality[57:], total_tonality[57:] > .02)),
        mean_specific_tonality_tu=_active_mean(spec_tonality[57:], active),
        tonal_frequency_hz=total_freq, specific_tonal_frequency_hz=freq,
        mean_specific_tonal_frequency_hz=_active_mean(freq[57:], active),
        tonal_loudness_sone_per_bark=tone, noise_loudness_sone_per_bark=noise,
        bark_axis=_BARK, centre_frequencies_hz=_CENTRES,
    )
    return EcmaTonalAnalysis(loudness=loudness, tonality=tonality)


def loudness_ecma(
    pressure_pa: ArrayLike, sample_rate_hz: int, *,
    sound_field: SoundField | str = SoundField.FREE,
) -> EcmaLoudnessResult:
    """Full Section 8 ECMA-418-2:2025 loudness in sone_HMS.

    Accepts calibrated mono or stereo pressure, at least 304 ms long.
    See ecma_tonal_analysis for resampling, timing and input conventions.
    Stereo results include the specific-pattern binaural quadratic mean.
    This differs from MoSQITo's Section 5 basis loudness.
    """
    return ecma_tonal_analysis(pressure_pa, sample_rate_hz, sound_field=sound_field).loudness


def tonality_ecma(
    pressure_pa: ArrayLike, sample_rate_hz: int, *,
    sound_field: SoundField | str = SoundField.FREE,
) -> EcmaTonalityResult:
    """Section 6 ECMA-418-2:2025 psychoacoustic tonality in tu_HMS.

    Accepts calibrated mono or stereo pressure, at least 304 ms long.
    See ecma_tonal_analysis for resampling, timing and input conventions.
    Stereo ears are evaluated independently; no channel mixing is applied.
    This is the Sottek model, distinct from Part 1 TNR/prominence ratio.
    """
    return ecma_tonal_analysis(pressure_pa, sample_rate_hz, sound_field=sound_field).tonality


def roughness_ecma(
    pressure_pa: ArrayLike, sample_rate_hz: int, *,
    sound_field: SoundField | str = SoundField.FREE,
) -> EcmaRoughnessResult:
    """Section 7 ECMA-418-2:2025 roughness in asper_HMS.

    Accepts calibrated mono (samples,) or stereo (samples, 2) pressure in
    pascals, at least 320 ms long. Audio rates from 8 to 192 kHz are resampled
    to 48 kHz. Output is at exactly 50 Hz, from zero through the last 20 ms
    grid point within the input duration. Representative roughness is the
    90th percentile after the first 320 ms. Stereo results include a
    binaural quadratic mean of the specific patterns. Optional entropy
    weighting (Section 7.1.6) is not applied. Inputs are never modified.
    """
    field = _field(sound_field)
    channels, mono = _signals(pressure_pa, sample_rate_hz, 15360)
    binding = get_native_library()
    frames = _frame_count(binding, "ms_ecma_roughness_frame_count",
                          binding.lib.ms_ecma_roughness_frame_count, channels[0].size)
    patterns = []
    for channel in channels:
        pattern = np.empty((frames, 53), dtype=np.float64)
        written = ctypes.c_size_t()
        status = binding.lib.ms_roughness_ecma(
            _pointer(channel), channel.size, NATIVE_SAMPLE_RATE_HZ, field,
            _pointer(pattern), pattern.size, ctypes.byref(written))
        _check(binding, "ms_roughness_ecma", status)
        patterns.append(pattern)
    specific = patterns[0] if mono else np.stack(patterns, axis=-1)
    total = .5 * specific.sum(axis=1)
    binaural = None if mono else np.hypot(specific[..., 0], specific[..., 1]) / np.sqrt(2)
    bin_total = None if binaural is None else .5 * binaural.sum(axis=1)
    return EcmaRoughnessResult(
        time_s=np.arange(frames, dtype=np.float64) / 50,
        roughness_asper=total, specific_roughness_asper_per_bark=specific,
        roughness_90_asper=_scalar_or_array(np.quantile(total[16:], .9, axis=0, method="hazen")),
        mean_specific_roughness_asper_per_bark=np.mean(specific[16:], axis=0),
        bark_axis=_BARK, centre_frequencies_hz=_CENTRES,
        binaural_roughness_asper=bin_total, binaural_specific_roughness_asper_per_bark=binaural,
        binaural_roughness_90_asper=None if bin_total is None else float(np.quantile(bin_total[16:], .9, method="hazen")),
    )
