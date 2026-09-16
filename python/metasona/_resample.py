# SPDX-License-Identifier: Apache-2.0
# MetaSona authors: Jiahua Zhang and Codex, September 2026.

"""The single sample-rate conversion boundary used by the Python package."""

from __future__ import annotations

import math

import numpy as np
from numpy.typing import NDArray
from scipy.signal import resample_poly

from .exceptions import MetaSonaValidationError

NATIVE_SAMPLE_RATE_HZ = 48_000
MAX_RATIONAL_FACTOR = 4_096


def to_native_rate(pressure_pa: NDArray[np.float64], sample_rate_hz: int) -> NDArray[np.float64]:
    """Convert mono pressure samples to exactly 48 kHz using polyphase FIR resampling.

    SciPy's polyphase resampler applies an anti-aliasing FIR filter. The rational
    conversion factors are reduced first, so common rates such as 44.1 kHz use
    the exact 160/147 ratio. Native kernels never receive any other sample rate.
    """
    if sample_rate_hz == NATIVE_SAMPLE_RATE_HZ:
        return np.array(pressure_pa, dtype=np.float64, copy=True, order="C")
    divisor = math.gcd(sample_rate_hz, NATIVE_SAMPLE_RATE_HZ)
    up = NATIVE_SAMPLE_RATE_HZ // divisor
    down = sample_rate_hz // divisor
    if max(up, down) > MAX_RATIONAL_FACTOR:
        raise MetaSonaValidationError(
            "sample_rate_hz produces an impractically large rational resampling filter; "
            "use a conventional audio rate or pre-resample explicitly"
        )
    try:
        converted = resample_poly(
            pressure_pa,
            up,
            down,
            window=("kaiser", 5.0),
            padtype="constant",
            cval=0.0,
        )
    except (MemoryError, OverflowError, ValueError) as exc:
        raise MetaSonaValidationError("polyphase resampling failed") from exc
    return np.ascontiguousarray(converted, dtype=np.float64)
