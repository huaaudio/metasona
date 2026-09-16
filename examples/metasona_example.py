# SPDX-FileCopyrightText: 2026 MetaSona contributors
# MetaSona author: Jiahua Zhang, September 2026.
# SPDX-License-Identifier: Apache-2.0

"""Calculate every public signal metric from calibrated mono pressure."""

from __future__ import annotations

import numpy as np

import metasona as ms


def calibrated_modulated_tone(sample_rate_hz: int, duration_s: float) -> np.ndarray:
    """Return a 1 kHz/70 Hz fully modulated signal at total 60 dB SPL."""
    count = round(sample_rate_hz * duration_s)
    time_s = np.arange(count, dtype=np.float64) / sample_rate_hz
    signal = (1.0 + np.sin(2.0 * np.pi * 70.0 * time_s)) * np.sin(2.0 * np.pi * 1_000.0 * time_s)
    target_rms_pa = 20e-6 * 10 ** (60.0 / 20.0)
    return signal * (target_rms_pa / np.sqrt(np.mean(np.square(signal))))


def main() -> None:
    # 44.1 kHz intentionally demonstrates the single Python resampling boundary.
    sample_rate_hz = 44_100
    pressure_pa = calibrated_modulated_tone(sample_rate_hz, 2.0)

    stationary = ms.stationary_loudness(
        pressure_pa,
        sample_rate_hz,
        sound_field=ms.SoundField.FREE,
    )
    varying = ms.time_varying_loudness(
        pressure_pa,
        sample_rate_hz,
        sound_field="free",
        time_skip_s=0.2,
    )
    roughness = ms.roughness_daniel_weber(
        pressure_pa,
        sample_rate_hz,
        time_skip_s=0.2,
    )
    tonality = ms.tonality_aures(
        pressure_pa,
        sample_rate_hz,
        sound_field="free",
        time_skip_s=0.2,
    )
    sharpness = ms.sharpness_din45692(stationary.specific_loudness_sone_per_bark)

    print(f"native library: {ms.native_version()}")
    print(f"stationary loudness: {stationary.loudness_sone:.6f} sone")
    print(f"sharpness: {sharpness.sharpness_acum:.6f} acum")
    print(
        f"time-varying loudness: {varying.time_s.size} frames, "
        f"mean {np.mean(varying.loudness_sone):.6f} sone"
    )
    print(
        f"roughness: {roughness.time_s.size} frames, "
        f"mean {np.mean(roughness.roughness_asper):.6f} asper"
    )
    print(f"tonality: {tonality.time_s.size} frames, mean {np.mean(tonality.tonality):.6f}")


if __name__ == "__main__":
    main()
