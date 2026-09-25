# SPDX-License-Identifier: Apache-2.0

from __future__ import annotations

import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import metasona as ms

from metasona import (
    LoudnessResult,
    MetaSonaValidationError,
    SharpnessResult,
)
from metasona._rolling import RollingAnalyzer, RollingMetric, RollingResult, RollingSnapshot


def _loudness_result(pressure_pa, sample_rate_hz, *, sound_field):
    del sample_rate_hz, sound_field
    value = float(np.mean(pressure_pa))
    return LoudnessResult(value, np.full(240, value), np.arange(240, dtype=np.float64))


def _snapshot_signature(snapshots):
    return [
        (
            snapshot.start_time_s,
            snapshot.end_time_s,
            snapshot.results[RollingMetric.STATIONARY_LOUDNESS].loudness_sone,
        )
        for snapshot in snapshots
    ]


class RollingPublicApiTests(unittest.TestCase):
    def test_rolling_names_are_exported_at_package_top_level(self):
        expected = {
            "RollingAnalyzer": RollingAnalyzer,
            "RollingMetric": RollingMetric,
            "RollingResult": RollingResult,
            "RollingSnapshot": RollingSnapshot,
        }

        for name, value in expected.items():
            with self.subTest(name=name):
                self.assertIs(getattr(ms, name), value)
                self.assertIn(name, ms.__all__)


class RollingConfigurationTests(unittest.TestCase):
    def test_accepts_enum_and_string_metrics(self):
        analyzer = RollingAnalyzer(
            48_000,
            [RollingMetric.STATIONARY_LOUDNESS, "tonality_aures"],
        )

        self.assertEqual(analyzer.push(np.zeros(1_000)), ())

    def test_rejects_invalid_metric_collections(self):
        invalid = (
            [],
            ["unknown"],
            [RollingMetric.STATIONARY_LOUDNESS, "stationary_loudness"],
            "stationary_loudness",
            None,
        )

        for metrics in invalid:
            with self.subTest(metrics=metrics):
                with self.assertRaises(MetaSonaValidationError):
                    RollingAnalyzer(48_000, metrics)

    def test_rejects_invalid_window_and_hop_values(self):
        for name, value in (
            ("window_s", 0.0),
            ("window_s", -1.0),
            ("window_s", float("nan")),
            ("window_s", float("inf")),
            ("window_s", True),
            ("window_s", "1"),
            ("hop_s", 0.0),
            ("hop_s", float("nan")),
        ):
            arguments = {name: value}
            with self.subTest(name=name, value=value):
                with self.assertRaises(MetaSonaValidationError):
                    RollingAnalyzer(
                        48_000,
                        [RollingMetric.STATIONARY_LOUDNESS],
                        **arguments,
                    )

        with self.assertRaises(MetaSonaValidationError):
            RollingAnalyzer(
                48_000,
                [RollingMetric.STATIONARY_LOUDNESS],
                window_s=0.5,
                hop_s=0.6,
            )

    def test_enforces_selected_metric_minimum_window(self):
        cases = (
            (RollingMetric.STATIONARY_LOUDNESS, 255 / 48_000),
            (RollingMetric.SHARPNESS_DIN45692, 255 / 48_000),
            (RollingMetric.ROUGHNESS_DANIEL_WEBER, 0.199),
            (RollingMetric.TONALITY_AURES, 0.249),
            (RollingMetric.LOUDNESS_ECMA, 0.303),
            (RollingMetric.TONALITY_ECMA, 0.303),
            (RollingMetric.ROUGHNESS_ECMA, 0.319),
        )

        for metric, window_s in cases:
            with self.subTest(metric=metric):
                with self.assertRaises(MetaSonaValidationError):
                    RollingAnalyzer(
                        48_000,
                        [metric],
                        window_s=window_s,
                        hop_s=window_s,
                    )

        RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.006,
            hop_s=0.006,
        )
        with self.assertRaises(MetaSonaValidationError):
            RollingAnalyzer(
                8_000,
                [RollingMetric.STATIONARY_LOUDNESS],
                window_s=0.005,
                hop_s=0.005,
            )

    def test_rejects_invalid_chunks_without_advancing_state(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.01,
        )

        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            self.assertEqual(analyzer.push(np.zeros(100)), ())
            for chunk in (
                np.array([], dtype=np.float64),
                np.array([np.nan]),
                np.zeros((2, 2)),
                np.array([True, False]),
                np.array([1 + 2j]),
            ):
                with self.subTest(chunk=chunk):
                    with self.assertRaises(MetaSonaValidationError):
                        analyzer.push(chunk)
            snapshots = analyzer.push(np.zeros(220))

        self.assertEqual(len(snapshots), 1)
        self.assertEqual(snapshots[0].end_time_s, 0.04)


class RollingSchedulingTests(unittest.TestCase):
    def test_emits_exact_windows_and_timestamps(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.01,
        )
        signal = np.arange(480, dtype=np.float64)

        with patch(
            "metasona._rolling.stationary_loudness",
            side_effect=_loudness_result,
        ) as stationary:
            self.assertEqual(analyzer.push(signal[:319]), ())
            snapshots = analyzer.push(signal[319:])

        self.assertEqual(
            [snapshot.start_time_s for snapshot in snapshots],
            [0.0, 0.01, 0.02],
        )
        self.assertEqual(
            [snapshot.end_time_s for snapshot in snapshots],
            [0.04, 0.05, 0.06],
        )
        expected_windows = (signal[:320], signal[80:400], signal[160:480])
        self.assertEqual(stationary.call_count, 3)
        for call, expected in zip(stationary.call_args_list, expected_windows):
            np.testing.assert_array_equal(call.args[0], expected)

    def test_chunk_partitioning_does_not_change_results(self):
        signal = np.arange(800, dtype=np.float64)

        def collect(parts):
            analyzer = RollingAnalyzer(
                8_000,
                [RollingMetric.STATIONARY_LOUDNESS],
                window_s=0.04,
                hop_s=0.02,
            )
            snapshots = []
            for part in parts:
                snapshots.extend(analyzer.push(part))
            return snapshots

        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            whole = collect([signal])
            partitioned = collect([signal[:123], signal[123:400], signal[400:]])

        self.assertEqual(_snapshot_signature(whole), _snapshot_signature(partitioned))

    def test_reset_discards_audio_and_restarts_timestamps(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.02,
        )

        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            first = analyzer.push(np.zeros(320))
            analyzer.reset()
            self.assertEqual(analyzer.push(np.zeros(319)), ())
            second = analyzer.push(np.zeros(1))

        self.assertEqual(first[0].start_time_s, 0.0)
        self.assertEqual(second[0].start_time_s, 0.0)
        self.assertEqual(second[0].end_time_s, 0.04)

    def test_retained_buffer_stays_below_one_window(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.01,
        )

        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            for _ in range(20):
                analyzer.push(np.zeros(100))
                self.assertLess(analyzer._buffer.size, 320)

    def test_snapshot_results_mapping_is_immutable(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.04,
        )

        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            snapshot = analyzer.push(np.zeros(320))[0]

        with self.assertRaises(TypeError):
            snapshot.results[RollingMetric.STATIONARY_LOUDNESS] = object()


class RollingDispatchTests(unittest.TestCase):
    def test_stationary_loudness_is_reused_for_sharpness(self):
        analyzer = RollingAnalyzer(
            48_000,
            [
                RollingMetric.STATIONARY_LOUDNESS,
                RollingMetric.SHARPNESS_DIN45692,
            ],
            window_s=0.01,
            hop_s=0.01,
        )
        loudness = _loudness_result(np.zeros(480), 48_000, sound_field="free")
        sharpness = SharpnessResult(1.25)

        with (
            patch("metasona._rolling.stationary_loudness", return_value=loudness) as stationary,
            patch("metasona._rolling.sharpness_din45692", return_value=sharpness) as sharp,
        ):
            snapshot = analyzer.push(np.zeros(480))[0]

        stationary.assert_called_once()
        sharp.assert_called_once_with(loudness.specific_loudness_sone_per_bark)
        self.assertIs(snapshot.results[RollingMetric.STATIONARY_LOUDNESS], loudness)
        self.assertIs(snapshot.results[RollingMetric.SHARPNESS_DIN45692], sharpness)

    def test_ecma_tonal_analysis_is_reused(self):
        analyzer = RollingAnalyzer(
            48_000,
            [RollingMetric.LOUDNESS_ECMA, RollingMetric.TONALITY_ECMA],
            window_s=0.304,
            hop_s=0.304,
        )
        loudness = object()
        tonality = object()
        combined = SimpleNamespace(loudness=loudness, tonality=tonality)

        with patch(
            "metasona._rolling.ecma_tonal_analysis",
            return_value=combined,
        ) as ecma:
            snapshot = analyzer.push(np.zeros(14_592))[0]

        ecma.assert_called_once()
        self.assertIs(snapshot.results[RollingMetric.LOUDNESS_ECMA], loudness)
        self.assertIs(snapshot.results[RollingMetric.TONALITY_ECMA], tonality)

    def test_independent_metric_dispatch(self):
        metrics = (
            RollingMetric.TIME_VARYING_LOUDNESS,
            RollingMetric.ROUGHNESS_DANIEL_WEBER,
            RollingMetric.TONALITY_AURES,
            RollingMetric.ROUGHNESS_ECMA,
        )
        analyzer = RollingAnalyzer(
            48_000,
            metrics,
            window_s=0.32,
            hop_s=0.32,
            sound_field="diffuse",
        )
        expected = {metric: object() for metric in metrics}

        with (
            patch(
                "metasona._rolling.time_varying_loudness",
                return_value=expected[RollingMetric.TIME_VARYING_LOUDNESS],
            ) as time_varying,
            patch(
                "metasona._rolling.roughness_daniel_weber",
                return_value=expected[RollingMetric.ROUGHNESS_DANIEL_WEBER],
            ) as daniel_weber,
            patch(
                "metasona._rolling.tonality_aures",
                return_value=expected[RollingMetric.TONALITY_AURES],
            ) as aures,
            patch(
                "metasona._rolling.roughness_ecma",
                return_value=expected[RollingMetric.ROUGHNESS_ECMA],
            ) as ecma_roughness,
        ):
            snapshot = analyzer.push(np.zeros(15_360))[0]

        time_varying.assert_called_once()
        daniel_weber.assert_called_once()
        aures.assert_called_once()
        ecma_roughness.assert_called_once()
        self.assertEqual(tuple(snapshot.results), metrics)
        for metric, result in expected.items():
            self.assertIs(snapshot.results[metric], result)

    def test_failed_computation_does_not_commit_state(self):
        analyzer = RollingAnalyzer(
            8_000,
            [RollingMetric.STATIONARY_LOUDNESS],
            window_s=0.04,
            hop_s=0.01,
        )
        signal = np.arange(400, dtype=np.float64)
        call_count = 0

        def fail_second_window(pressure_pa, sample_rate_hz, *, sound_field):
            nonlocal call_count
            call_count += 1
            if call_count == 2:
                raise RuntimeError("calculation failed")
            return _loudness_result(pressure_pa, sample_rate_hz, sound_field=sound_field)

        with patch(
            "metasona._rolling.stationary_loudness",
            side_effect=fail_second_window,
        ):
            with self.assertRaisesRegex(RuntimeError, "calculation failed"):
                analyzer.push(signal)

        self.assertEqual(analyzer._samples_seen, 0)
        self.assertEqual(analyzer._buffer.size, 0)
        with patch("metasona._rolling.stationary_loudness", _loudness_result):
            snapshots = analyzer.push(signal)
        self.assertEqual([item.end_time_s for item in snapshots], [0.04, 0.05])


if __name__ == "__main__":
    unittest.main()
