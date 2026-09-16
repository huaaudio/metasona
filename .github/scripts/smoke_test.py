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

    stationary = ms.stationary_loudness(pressure, sample_rate)
    assert np.isfinite(stationary.loudness_sone) and stationary.loudness_sone > 0
    specific = stationary.specific_loudness_sone_per_bark
    assert np.all(np.isfinite(specific)) and not specific.flags.writeable
    sharpness = ms.sharpness_din45692(specific).sharpness_acum
    assert np.isfinite(sharpness) and sharpness > 0

    varying = ms.time_varying_loudness(pressure, sample_rate)
    roughness = ms.roughness_daniel_weber(pressure, sample_rate)
    tonality = ms.tonality_aures(pressure, sample_rate)
    for timestamps, values in (
        (varying.time_s, varying.loudness_sone),
        (roughness.time_s, roughness.roughness_asper),
        (tonality.time_s, tonality.tonality),
    ):
        assert values.size > 0 and values.shape == timestamps.shape
        assert np.all(np.isfinite(values)) and np.all(values >= 0)
        assert np.all(np.isfinite(timestamps)) and np.all(np.diff(timestamps) > 0)
        assert not values.flags.writeable
    print(f"MetaSona {version('metasona')}: installed native library smoke test passed")


if __name__ == "__main__":
    main()
