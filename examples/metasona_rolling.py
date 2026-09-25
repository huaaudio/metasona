# SPDX-License-Identifier: Apache-2.0

"""Feed calibrated audio chunks through the rolling analysis API."""

from __future__ import annotations

import numpy as np

import metasona as ms


def calibrated_modulated_tone(sample_rate_hz: int, duration_s: float) -> np.ndarray:
    """Return a 1 kHz/70 Hz fully modulated signal at total 60 dB SPL."""
    sample_count = round(sample_rate_hz * duration_s)
    time_s = np.arange(sample_count, dtype=np.float64) / sample_rate_hz
    signal = (1.0 + np.sin(2.0 * np.pi * 70.0 * time_s)) * np.sin(
        2.0 * np.pi * 1_000.0 * time_s
    )
    target_rms_pa = 20e-6 * 10 ** (60.0 / 20.0)
    return signal * (target_rms_pa / np.sqrt(np.mean(np.square(signal))))


def main() -> None:
    sample_rate_hz = 48_000
    pressure_pa = calibrated_modulated_tone(sample_rate_hz, 2.0)
    analyzer = ms.RollingAnalyzer(
        sample_rate_hz,
        [
            ms.RollingMetric.STATIONARY_LOUDNESS,
            ms.RollingMetric.TONALITY_AURES,
        ],
        window_s=0.5,
        hop_s=0.25,
    )

    chunk_sizes = (4_096, 7_000, 2_048)
    offset = 0
    chunk_index = 0
    while offset < pressure_pa.size:
        chunk_size = chunk_sizes[chunk_index % len(chunk_sizes)]
        chunk = pressure_pa[offset : offset + chunk_size]
        offset += chunk.size
        chunk_index += 1
        for snapshot in analyzer.push(chunk):
            loudness = snapshot.results[ms.RollingMetric.STATIONARY_LOUDNESS]
            tonality = snapshot.results[ms.RollingMetric.TONALITY_AURES]
            print(
                f"{snapshot.end_time_s:5.2f} s  "
                f"{loudness.loudness_sone:7.3f} sone  "
                f"{np.mean(tonality.tonality):7.3f} mean tonality"
            )


if __name__ == "__main__":
    main()