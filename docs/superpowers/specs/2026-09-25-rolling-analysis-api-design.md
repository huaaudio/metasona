# Rolling Analysis API Design

## Purpose

MetaSona 0.2.2 will provide a reusable Python API for applications that receive
calibrated audio incrementally and want periodic sound-quality calculations.
MetaSona remains a computation library: it will not capture audio, create
threads, own an event loop, or render a user interface.

The API will adapt the existing batch calculations through bounded rolling
windows. It will accurately describe this behavior as rolling analysis rather
than claiming that the native kernels preserve state between calls.

## Goals

- Accept arbitrary-sized chunks of calibrated mono pressure samples.
- Let callers select any supported signal metric.
- Emit immutable, timestamped results at a sample-accurate cadence.
- Return the existing typed result objects without inventing new scalar
  definitions for live display.
- Reuse shared calculations within one snapshot.
- Add no runtime, native, audio-capture, or UI dependencies.
- Preserve Python 3.9 through 3.14 support and native ABI version 1.

## Non-goals

- Audio-device discovery or microphone capture.
- Background workers, callbacks, asyncio integration, or scheduling policy.
- Hard real-time guarantees or use from an audio callback thread.
- Incremental native filter state across windows.
- Stereo rolling analysis in 0.2.2. Existing ECMA batch APIs retain stereo
  support.
- Dropping work when computation falls behind.
- Padding or emitting incomplete final windows.

## Public API

Four names will be exported from `metasona`:

```python
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
    start_time_s: float
    end_time_s: float
    results: Mapping[RollingMetric, RollingResult]


class RollingAnalyzer:
    def __init__(
        self,
        sample_rate_hz: int,
        metrics: Iterable[RollingMetric | str],
        *,
        window_s: float = 1.0,
        hop_s: float = 0.2,
        sound_field: SoundField | str = SoundField.FREE,
    ) -> None: ...

    def push(self, pressure_pa: ArrayLike) -> tuple[RollingSnapshot, ...]: ...
    def reset(self) -> None: ...
```

`RollingSnapshot.results` is an immutable mapping. Its values remain the
existing frozen result objects with read-only arrays.

Typical use:

```python
analyzer = metasona.RollingAnalyzer(
    sample_rate_hz=48_000,
    metrics=[
        metasona.RollingMetric.STATIONARY_LOUDNESS,
        metasona.RollingMetric.TONALITY_AURES,
    ],
    window_s=1.0,
    hop_s=0.2,
)

for snapshot in analyzer.push(pressure_chunk_pa):
    loudness = snapshot.results[
        metasona.RollingMetric.STATIONARY_LOUDNESS
    ]
```

## Input and Configuration

The sample rate is fixed for an analyzer instance and uses the existing
integer range of 8 through 192 kHz. Each pushed chunk must be a finite,
nonempty, one-dimensional real numeric array. The analyzer copies input data;
callers may safely reuse or mutate their buffers after `push` returns.

`metrics` must be a nonempty iterable of enum values or their exact string
values. Unknown and duplicate metrics are configuration errors. Metric order
does not affect execution or output mapping semantics.

`window_s` and `hop_s` must be finite positive real numbers, excluding
booleans. They convert once at construction to integer sample counts with
`round(value * sample_rate_hz)`. Both rounded counts must be positive and the
hop must not exceed the window. The rounded sample counts control all
scheduling; the original floating-point arguments are not used afterward.

The rolling window must meet the largest minimum selected by its metrics:

| Metric | Minimum native window |
|---|---:|
| Stationary loudness / sharpness | 256 samples at 48 kHz |
| Time-varying loudness | 1 sample at 48 kHz |
| Daniel-Weber roughness | 200 ms |
| Aures tonality | 250 ms |
| ECMA loudness / tonality | 304 ms |
| ECMA roughness | 320 ms |

Minimum validation uses integer sample-rate arithmetic so rounding and
resampling cannot produce a window shorter than the native requirement.
`sound_field` is validated once and applies only to metrics that already use
it; selecting another metric does not make it an error.

## Scheduling and Buffering

The analyzer maintains an absolute count of accepted input samples. The first
snapshot ends at exactly `window_samples`; subsequent snapshots end every
`hop_samples`. A call to `push` returns an empty tuple when no window has
completed and may return multiple snapshots when a large chunk crosses
multiple boundaries.

For a snapshot ending at absolute sample `end`, the analyzer passes exactly
`[end - window_samples:end]` to each selected batch calculation. Snapshot
start and end timestamps are these absolute sample positions divided by the
configured sample rate. Time axes nested inside existing result objects remain
relative to the rolling window, as they are in the batch APIs.

After successful processing, samples that cannot participate in the next
window are discarded. Retained internal audio therefore remains below one
window except while the current input chunk is being processed. There is no
unbounded result queue.

`reset()` discards buffered samples and restarts the sample clock and schedule
at zero. It emits no partial result. Applications must call it after stream
discontinuities or sample loss.

## Metric Dispatch and Reuse

Each selected metric maps to the corresponding existing public calculation:

- `STATIONARY_LOUDNESS` -> `stationary_loudness`
- `TIME_VARYING_LOUDNESS` -> `time_varying_loudness`
- `ROUGHNESS_DANIEL_WEBER` -> `roughness_daniel_weber`
- `TONALITY_AURES` -> `tonality_aures`
- `SHARPNESS_DIN45692` -> `sharpness_din45692`
- `LOUDNESS_ECMA` -> `loudness_ecma`
- `TONALITY_ECMA` -> `tonality_ecma`
- `ROUGHNESS_ECMA` -> `roughness_ecma`

Within a snapshot, sharpness reuses the specific-loudness pattern from one
stationary-loudness calculation. If stationary loudness was not selected, its
intermediate result is not exposed. Likewise, selecting either or both ECMA
loudness and tonality invokes `ecma_tonal_analysis` once; only selected result
objects are inserted into the mapping.

All calculations are synchronous. The analyzer computes every due snapshot in
order and never silently skips one. Applications with expensive selections
must invoke `push` outside the real-time capture callback, for example from
their own worker thread.

## Errors and State

Configuration and chunk validation failures raise `MetaSonaValidationError`.
Native failures retain their existing `NativeCallError` and
`NativeLibraryError` behavior.

`push` is transactional with respect to analyzer state. It validates the whole
chunk and computes all due snapshots using temporary scheduling state. If any
calculation raises, no input, sample count, or next-boundary update is
committed, and no partial snapshot tuple is returned. The caller may retry the
same chunk after handling the failure.

Instances are not thread-safe. One application-owned execution context must
serialize calls to `push` and `reset`.

## Performance Expectations

A September 2026 Windows benchmark on a one-second, 48 kHz mono window measured
approximately:

| Calculation | Median time |
|---|---:|
| Stationary loudness | 9.5 ms |
| Time-varying loudness | 21.2 ms |
| Aures tonality | 41.4 ms |
| ECMA roughness | 138.3 ms |
| Daniel-Weber roughness | 158.0 ms |
| Combined ECMA loudness/tonality | 615.2 ms |

These are workload observations, not timing guarantees. Documentation will
explain that callers choose a hop compatible with their selected metrics and
hardware. The API will not reject an aggressive hop because applications may
process faster hardware, tolerate backlog, or submit chunks from stored data.

One second of mono `float64` audio at 48 kHz occupies about 384 KiB. At the
maximum accepted 192 kHz it occupies about 1.5 MiB. Existing calculations may
allocate larger temporary workspaces independently of the rolling buffer.

## Dependencies and Compatibility

The implementation is Python-only and uses NumPy already required by
MetaSona. Existing batch functions continue using SciPy for resampling. No new
runtime or build dependency is introduced. In particular, MetaSona will not
depend on audio-device, terminal UI, web, async, or concurrency packages.

No C header, exported symbol, native implementation, ABI version, or wheel
platform changes. Existing public APIs retain their signatures and behavior.

## Files and Responsibilities

- `python/metasona/_rolling.py`: metric enum, result alias, immutable snapshot,
  configuration validation, buffering, scheduling, and metric dispatch.
- `python/metasona/__init__.py`: public exports.
- `tests/test_rolling.py`: scheduling, validation, chunk partitioning,
  transactional behavior, calculation reuse, reset, and immutability tests
  using `unittest` and existing dependencies.
- `.github/scripts/smoke_test.py`: installed-wheel end-to-end rolling check.
- `docs/python-api.md`: public contract and performance guidance.
- `examples/metasona_rolling.py`: device-independent synthetic chunk example.
- `README.md`: concise feature mention and link to API documentation.

## Verification

Behavioral tests will cover:

- no output before the first complete window;
- exact timestamps and multiple results from one large chunk;
- identical snapshots regardless of input chunk partitioning;
- selected result keys and existing result types;
- reuse of stationary and ECMA combined calculations;
- finite/shape/configuration validation without state mutation;
- computation failure without state mutation;
- bounded retained audio after repeated pushes;
- reset semantics;
- immutable snapshot mapping and result arrays;
- Python 3.9-compatible syntax and imports;
- an installed wheel loading its bundled native library and producing rolling
  results through the existing cross-platform publish smoke test.

## Estimated Effort

| Area | Estimated changed LOC |
|---|---:|
| Rolling API and validation | 260-340 |
| Exports and integration | 15-30 |
| Behavioral tests | 280-420 |
| Installed-package smoke coverage | 20-40 |
| Documentation and example | 110-160 |
| **Total** | **685-990** |

No native C lines are expected. Supporting stereo conditionally by metric
would add roughly 100-180 lines of implementation and tests and is deferred.
A true stateful native streaming ABI would instead require approximately
1,800-3,000 lines across native state objects, lifecycle functions, wrappers,
tests, and documentation.

## Alternatives Rejected

### Background worker and callbacks

This would move scheduling policy into a base computation library and add
queue capacity, cancellation, shutdown, exception delivery, and thread-safety
contracts. Applications already have their own execution models, so this
complexity belongs above MetaSona.

### True incremental native kernels

This could reduce overlapping recomputation and preserve continuous filter
state, but every metric has a different state and framing model. It would
expand the C ABI and require substantial reference-equivalence work. It is not
appropriate for a dependency-fix patch release.

### Scalar-only snapshots

Choosing one scalar per metric would add new semantic decisions, especially
for time-varying and ECMA outputs. Returning existing result objects keeps the
rolling API compositional and lets each application select its own display or
aggregation statistic.
