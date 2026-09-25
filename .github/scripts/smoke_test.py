"""Exercise the installed distribution and its bundled native library."""

from importlib.metadata import version
from pathlib import Path

import numpy as np

import metasona as ms


def main() -> None:
    checkout = Path(__file__).resolve().parents[2]
    assert not Path(ms.__file__).resolve().is_relative_to(checkout / "python")
    assert ms.native_version() == version("metasona")

    # Exercise resampling as well as every public signal metric.
    sample_rate = 44_100
    time = np.arange(sample_rate, dtype=np.float64) / sample_rate
    pressure = (1 + np.sin(2 * np.pi * 70 * time)) * np.sin(2 * np.pi * 1000 * time)
    pressure *= 0.02 / np.sqrt(np.mean(pressure**2))

    rolling_metrics = {
        ms.RollingMetric.STATIONARY_LOUDNESS,
        ms.RollingMetric.SHARPNESS_DIN45692,
        ms.RollingMetric.LOUDNESS_ECMA,
        ms.RollingMetric.TONALITY_ECMA,
    }
    rolling = ms.RollingAnalyzer(
        sample_rate,
        rolling_metrics,
        window_s=0.5,
        hop_s=0.5,
    )
    snapshots = []
    for chunk in (pressure[:12_345], pressure[12_345:30_000], pressure[30_000:]):
        snapshots.extend(rolling.push(chunk))
    assert [snapshot.end_time_s for snapshot in snapshots] == [0.5, 1.0]
    for snapshot in snapshots:
        assert set(snapshot.results) == rolling_metrics
        assert np.isfinite(snapshot.results[ms.RollingMetric.STATIONARY_LOUDNESS].loudness_sone)
        assert np.isfinite(snapshot.results[ms.RollingMetric.SHARPNESS_DIN45692].sharpness_acum)
        assert np.isfinite(snapshot.results[ms.RollingMetric.LOUDNESS_ECMA].mean_loudness_sone)
        assert np.isfinite(snapshot.results[ms.RollingMetric.TONALITY_ECMA].mean_tonality_tu)

    stationary = ms.stationary_loudness(pressure, sample_rate)
    assert np.isfinite(stationary.loudness_sone) and stationary.loudness_sone > 0
    specific = stationary.specific_loudness_sone_per_bark
    assert np.all(np.isfinite(specific)) and not specific.flags.writeable
    sharpness = ms.sharpness_din45692(specific).sharpness_acum
    assert np.isfinite(sharpness) and sharpness > 0

    varying = ms.time_varying_loudness(pressure, sample_rate)
    roughness = ms.roughness_daniel_weber(pressure, sample_rate)
    tonality = ms.tonality_aures(pressure, sample_rate)
    ecma = ms.ecma_tonal_analysis(pressure, sample_rate)
    ecma_roughness = ms.roughness_ecma(pressure, sample_rate)
    assert 0.97 < ecma_roughness.roughness_90_asper < 1.03
    tone = np.sqrt(2) * .002 * np.sin(2 * np.pi * 1000 * time)
    calibration = ms.ecma_tonal_analysis(tone, sample_rate)
    assert abs(calibration.loudness.mean_loudness_sone - 1) < .01
    assert abs(calibration.tonality.mean_tonality_tu - 1) < .01
    for timestamps, values in (
        (varying.time_s, varying.loudness_sone),
        (roughness.time_s, roughness.roughness_asper),
        (tonality.time_s, tonality.tonality),
        (ecma.loudness.time_s, ecma.loudness.loudness_sone),
        (ecma.tonality.time_s, ecma.tonality.tonality_tu),
        (ecma_roughness.time_s, ecma_roughness.roughness_asper),
    ):
        assert values.size > 0 and values.shape == timestamps.shape
        assert np.all(np.isfinite(values)) and np.all(values >= 0)
        assert np.all(np.isfinite(timestamps)) and np.all(np.diff(timestamps) > 0)
        assert not values.flags.writeable
    print(f"MetaSona {version('metasona')}: installed native library smoke test passed")


if __name__ == "__main__":
    main()
