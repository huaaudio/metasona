# SPDX-License-Identifier: Apache-2.0

"""Synchronous rolling-window analysis over incrementally supplied audio."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
from numbers import Real
from types import MappingProxyType
from typing import Union

import numpy as np
from numpy.typing import ArrayLike

from ._api import (
    roughness_daniel_weber,
    sharpness_din45692,
    stationary_loudness,
    time_varying_loudness,
    tonality_aures,
)
from ._ecma import ecma_tonal_analysis, roughness_ecma
from ._ecma_types import EcmaLoudnessResult, EcmaRoughnessResult, EcmaTonalityResult
from ._types import (
    LoudnessResult,
    RoughnessResult,
    SharpnessResult,
    SoundField,
    TimeVaryingLoudnessResult,
    TonalityResult,
)
from ._validation import mono_pressure, positive_sample_rate, sound_field as parse_sound_field
from .exceptions import MetaSonaValidationError

RollingResult = Union[
    LoudnessResult,
    TimeVaryingLoudnessResult,
    RoughnessResult,
    TonalityResult,
    SharpnessResult,
    EcmaLoudnessResult,
    EcmaTonalityResult,
    EcmaRoughnessResult,
]


class RollingMetric(str, Enum):
    """Signal metric available through :class:`RollingAnalyzer`."""

    STATIONARY_LOUDNESS = "stationary_loudness"
    TIME_VARYING_LOUDNESS = "time_varying_loudness"
    ROUGHNESS_DANIEL_WEBER = "roughness_daniel_weber"
    TONALITY_AURES = "tonality_aures"
    SHARPNESS_DIN45692 = "sharpness_din45692"
    LOUDNESS_ECMA = "loudness_ecma"
    TONALITY_ECMA = "tonality_ecma"
    ROUGHNESS_ECMA = "roughness_ecma"


@dataclass(frozen=True)
class RollingSnapshot:
    """Results for one complete rolling window on the stream timeline."""

    start_time_s: float
    end_time_s: float
    results: Mapping[RollingMetric, RollingResult]

    def __post_init__(self) -> None:
        object.__setattr__(self, "results", MappingProxyType(dict(self.results)))


_MINIMUM_NATIVE_SAMPLES = {
    RollingMetric.STATIONARY_LOUDNESS: 256,
    RollingMetric.TIME_VARYING_LOUDNESS: 1,
    RollingMetric.ROUGHNESS_DANIEL_WEBER: 9_600,
    RollingMetric.TONALITY_AURES: 12_000,
    RollingMetric.SHARPNESS_DIN45692: 256,
    RollingMetric.LOUDNESS_ECMA: 14_592,
    RollingMetric.TONALITY_ECMA: 14_592,
    RollingMetric.ROUGHNESS_ECMA: 15_360,
}


def _parse_metrics(values: Iterable[RollingMetric | str]) -> tuple[RollingMetric, ...]:
    if values is None or isinstance(values, (str, RollingMetric)):
        raise MetaSonaValidationError("metrics must be a nonempty iterable of metric values")
    try:
        supplied = list(values)
    except TypeError as exc:
        raise MetaSonaValidationError(
            "metrics must be a nonempty iterable of metric values"
        ) from exc
    if not supplied:
        raise MetaSonaValidationError("metrics must not be empty")

    parsed = []
    seen = set()
    for value in supplied:
        if isinstance(value, RollingMetric):
            metric = value
        elif isinstance(value, str):
            try:
                metric = RollingMetric(value)
            except ValueError as exc:
                raise MetaSonaValidationError(f"unknown rolling metric: {value!r}") from exc
        else:
            raise MetaSonaValidationError(f"invalid rolling metric: {value!r}")
        if metric in seen:
            raise MetaSonaValidationError(f"duplicate rolling metric: {metric.value!r}")
        parsed.append(metric)
        seen.add(metric)
    return tuple(parsed)


def _duration_samples(value: Real, sample_rate_hz: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise MetaSonaValidationError(f"{name} must be a finite positive number")
    duration_s = float(value)
    if not math.isfinite(duration_s) or duration_s <= 0.0:
        raise MetaSonaValidationError(f"{name} must be a finite positive number")
    samples = round(duration_s * sample_rate_hz)
    if samples <= 0:
        raise MetaSonaValidationError(f"{name} rounds to zero samples")
    return samples


class RollingAnalyzer:
    """Compute selected MetaSona metrics over complete rolling windows.

    The analyzer performs no background work and is not thread-safe. Callers
    own capture, scheduling, and serialization of ``push`` and ``reset``.
    """

    def __init__(
        self,
        sample_rate_hz: int,
        metrics: Iterable[RollingMetric | str],
        *,
        window_s: float = 1.0,
        hop_s: float = 0.2,
        sound_field: SoundField | str = SoundField.FREE,
    ) -> None:
        self._sample_rate_hz = positive_sample_rate(sample_rate_hz)
        self._metrics = _parse_metrics(metrics)
        self._window_samples = _duration_samples(
            window_s, self._sample_rate_hz, "window_s"
        )
        self._hop_samples = _duration_samples(hop_s, self._sample_rate_hz, "hop_s")
        if self._hop_samples > self._window_samples:
            raise MetaSonaValidationError("hop_s must not exceed window_s")

        minimum_native = max(_MINIMUM_NATIVE_SAMPLES[metric] for metric in self._metrics)
        if self._window_samples * 48_000 < minimum_native * self._sample_rate_hz:
            minimum_s = minimum_native / 48_000
            raise MetaSonaValidationError(
                f"window_s is too short for the selected metrics; use at least {minimum_s:g} s"
            )
        self._sound_field = parse_sound_field(sound_field)
        self.reset()

    def reset(self) -> None:
        """Discard buffered audio and restart the stream timeline at zero."""
        self._buffer = np.empty(0, dtype=np.float64)
        self._buffer_start = 0
        self._samples_seen = 0
        self._next_end = self._window_samples

    def push(self, pressure_pa: ArrayLike) -> tuple[RollingSnapshot, ...]:
        """Accept one mono chunk and return every newly completed snapshot."""
        chunk = mono_pressure(pressure_pa)
        combined = np.concatenate((self._buffer, chunk))
        combined_end = self._samples_seen + chunk.size
        next_end = self._next_end
        snapshots = []

        while next_end <= combined_end:
            window_start = next_end - self._window_samples
            relative_start = window_start - self._buffer_start
            relative_end = next_end - self._buffer_start
            window = combined[relative_start:relative_end]
            results = self._calculate(window)
            snapshots.append(
                RollingSnapshot(
                    window_start / self._sample_rate_hz,
                    next_end / self._sample_rate_hz,
                    results,
                )
            )
            next_end += self._hop_samples

        retained_start = next_end - self._window_samples
        relative_retained = retained_start - self._buffer_start
        retained = np.array(combined[relative_retained:], dtype=np.float64, copy=True, order="C")

        self._buffer = retained
        self._buffer_start = retained_start
        self._samples_seen = combined_end
        self._next_end = next_end
        return tuple(snapshots)

    def _calculate(self, window: np.ndarray) -> Mapping[RollingMetric, RollingResult]:
        available = {}
        metric_set = set(self._metrics)

        loudness_metrics = {
            RollingMetric.STATIONARY_LOUDNESS,
            RollingMetric.SHARPNESS_DIN45692,
        }
        if metric_set & loudness_metrics:
            loudness = stationary_loudness(
                window,
                self._sample_rate_hz,
                sound_field=self._sound_field,
            )
            if RollingMetric.STATIONARY_LOUDNESS in metric_set:
                available[RollingMetric.STATIONARY_LOUDNESS] = loudness
            if RollingMetric.SHARPNESS_DIN45692 in metric_set:
                available[RollingMetric.SHARPNESS_DIN45692] = sharpness_din45692(
                    loudness.specific_loudness_sone_per_bark
                )

        if RollingMetric.TIME_VARYING_LOUDNESS in metric_set:
            available[RollingMetric.TIME_VARYING_LOUDNESS] = time_varying_loudness(
                window,
                self._sample_rate_hz,
                sound_field=self._sound_field,
            )
        if RollingMetric.ROUGHNESS_DANIEL_WEBER in metric_set:
            available[RollingMetric.ROUGHNESS_DANIEL_WEBER] = roughness_daniel_weber(
                window,
                self._sample_rate_hz,
            )
        if RollingMetric.TONALITY_AURES in metric_set:
            available[RollingMetric.TONALITY_AURES] = tonality_aures(
                window,
                self._sample_rate_hz,
                sound_field=self._sound_field,
            )

        ecma_metrics = {RollingMetric.LOUDNESS_ECMA, RollingMetric.TONALITY_ECMA}
        if metric_set & ecma_metrics:
            analysis = ecma_tonal_analysis(
                window,
                self._sample_rate_hz,
                sound_field=self._sound_field,
            )
            if RollingMetric.LOUDNESS_ECMA in metric_set:
                available[RollingMetric.LOUDNESS_ECMA] = analysis.loudness
            if RollingMetric.TONALITY_ECMA in metric_set:
                available[RollingMetric.TONALITY_ECMA] = analysis.tonality
        if RollingMetric.ROUGHNESS_ECMA in metric_set:
            available[RollingMetric.ROUGHNESS_ECMA] = roughness_ecma(
                window,
                self._sample_rate_hz,
                sound_field=self._sound_field,
            )

        return {metric: available[metric] for metric in self._metrics}