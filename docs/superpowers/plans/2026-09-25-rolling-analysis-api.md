# Rolling Analysis API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a synchronous, selectable-metric rolling analyzer that turns incrementally supplied calibrated mono audio into timestamped existing MetaSona result objects.

**Architecture:** A Python-only `RollingAnalyzer` owns a bounded input-rate buffer and schedules exact sample-count rolling windows. It dispatches each completed window through existing batch APIs, reuses stationary and ECMA combined calculations, and commits state only after all due snapshots succeed.

**Tech Stack:** Python 3.9-3.14, NumPy, existing MetaSona Python/native APIs, `unittest`

**Spec:** `docs/superpowers/specs/2026-09-25-rolling-analysis-api-design.md`

## Global Constraints

- Add no runtime, build, capture, UI, threading, or async dependencies.
- Preserve Python 3.9 through 3.14 support.
- Preserve native ABI version 1 and make no native C changes.
- Accept calibrated mono pressure only; existing ECMA batch stereo support is unchanged.
- All scheduling uses rounded integer sample counts and emits no incomplete window.
- `push` is synchronous and transactional; instances are not thread-safe.

## Review Focus

- Chunk boundaries must not alter windows, timestamps, or numerical results.
- Failed validation or computation must not advance accepted samples or scheduling state.
- A large push may produce several ordered snapshots without retaining obsolete audio.
- Sharpness and paired ECMA selections must not duplicate their shared calculation.
- Public mappings and nested result arrays must remain immutable.

---

### Task 1: Core Rolling Analyzer

**Files:**
- Create: `tests/test_rolling.py`
- Create: `python/metasona/_rolling.py`

**Interfaces:**
- Consumes: existing metric functions, result classes, `mono_pressure`, `positive_sample_rate`, `sound_field`, and `MetaSonaValidationError`
- Produces: `RollingResult`, `RollingMetric`, `RollingSnapshot`, and `RollingAnalyzer.push/reset`

- [ ] **Step 1: Write failing scheduling and configuration tests**

Create `tests/test_rolling.py` with `unittest` cases that import the four new names and verify:

```python
analyzer = RollingAnalyzer(
    8_000,
    [RollingMetric.STATIONARY_LOUDNESS],
    window_s=0.04,
    hop_s=0.01,
)
assert analyzer.push(np.zeros(319)) == ()
snapshots = analyzer.push(np.zeros(161))
assert [item.end_time_s for item in snapshots] == [0.04, 0.05, 0.06]
```

Patch `_rolling.stationary_loudness` with a deterministic real
`LoudnessResult` factory. Add cases for unknown/duplicate/empty metrics,
invalid duration values, `hop_s > window_s`, metric-specific minimum windows,
and invalid chunks.

- [ ] **Step 2: Run tests to verify RED**

Run:

```powershell
$env:PYTHONPATH="$PWD/python"
python -m unittest discover -s tests -p test_rolling.py -v
```

Expected: FAIL because `metasona._rolling` does not exist.

- [ ] **Step 3: Implement public rolling types and validation**

Create `_rolling.py` with:

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

class RollingMetric(str, Enum): ...

@dataclass(frozen=True)
class RollingSnapshot:
    start_time_s: float
    end_time_s: float
    results: Mapping[RollingMetric, RollingResult]

    def __post_init__(self) -> None:
        object.__setattr__(self, "results", MappingProxyType(dict(self.results)))
```

Implement strict metric parsing, finite positive duration conversion, hop/window checks, metric minimum checks using
`window_samples * 48_000 >= native_minimum * sample_rate_hz`, and initialized empty state.

- [ ] **Step 4: Implement transactional scheduling and dispatch**

Implement `push` by validating/copying the chunk, concatenating with retained
audio in local state, computing each due absolute window, and only then
committing `_buffer`, `_buffer_start`, `_samples_seen`, and `_next_end`.
Build selected output in configured metric order. Compute stationary loudness
once when loudness or sharpness is selected, and call `ecma_tonal_analysis`
once when ECMA loudness or tonality is selected. Implement `reset` by restoring
constructor scheduling state.

- [ ] **Step 5: Expand tests for complete behavior**

Add cases that verify:

```python
# One chunk and uneven chunks produce the same snapshot times and values.
whole = collect([signal])
partitioned = collect([signal[:123], signal[123:400], signal[400:]])
self.assertEqual(snapshot_signature(whole), snapshot_signature(partitioned))
```

Also verify multiple due windows, reset, mapping immutability, bounded retained
buffer, exact selected keys/result objects, one stationary call for
loudness+sharpness, one ECMA call for loudness+tonality, and retrying the same
chunk after a patched calculation fails midway.

- [ ] **Step 6: Run focused tests to verify GREEN**

Run the Task 1 command again.
Expected: all rolling tests pass.

- [ ] **Step 7: Commit Task 1**

```powershell
git add python/metasona/_rolling.py tests/test_rolling.py
git commit -m "Add rolling sound metric analyzer"
```

### Task 2: Public Package and Distribution Integration

**Files:**
- Modify: `python/metasona/__init__.py`
- Modify: `.github/scripts/smoke_test.py`
- Test: `tests/test_rolling.py`

**Interfaces:**
- Consumes: the four names from Task 1
- Produces: stable top-level imports and installed-wheel smoke coverage

- [ ] **Step 1: Write failing top-level export test**

Add a test importing `RollingAnalyzer`, `RollingMetric`, `RollingResult`, and
`RollingSnapshot` directly from `metasona`, and assert every name appears in
`metasona.__all__`.

- [ ] **Step 2: Run test to verify RED**

Run the Task 1 test command.
Expected: FAIL because top-level exports are absent.

- [ ] **Step 3: Export the public API**

Import the four names from `._rolling` in `__init__.py` and add them to
`__all__`, preserving the existing import-time guarantee that native code is
not loaded until a calculation runs.

- [ ] **Step 4: Extend installed-distribution smoke coverage**

In `.github/scripts/smoke_test.py`, feed its existing calibrated pressure in
uneven chunks to a rolling analyzer selecting stationary loudness, sharpness,
ECMA loudness, and ECMA tonality. Assert snapshots are produced, selected keys
are exact, timestamps increase, and values are finite. This exercises the
bundled native library and shared-calculation paths in sdist and every wheel.

- [ ] **Step 5: Run integration checks**

Run:

```powershell
$env:PYTHONPATH="$PWD/python"
$env:METASONA_LIBRARY="$PWD/build/py3-none-win_amd64/metasona.dll"
python -m unittest discover -s tests -p test_rolling.py -v
python .github/scripts/smoke_test.py
```

Expected: all unit tests and installed/source smoke assertions pass. If the
smoke script's installed-package guard rejects source execution, build and
install the wheel into a temporary uv environment before rerunning it.

- [ ] **Step 6: Commit Task 2**

```powershell
git add python/metasona/__init__.py .github/scripts/smoke_test.py tests/test_rolling.py
git commit -m "Expose rolling analyzer API"
```

### Task 3: Documentation, Example, and Release Verification

**Files:**
- Create: `examples/metasona_rolling.py`
- Modify: `docs/python-api.md`
- Modify: `README.md`

**Interfaces:**
- Consumes: top-level API from Task 2
- Produces: device-independent usage guidance and copy-paste example

- [ ] **Step 1: Add the synthetic chunk example**

Create an example that generates two seconds of calibrated modulated mono audio,
constructs a `RollingAnalyzer` with stationary loudness and Aures tonality,
pushes uneven chunks, and prints snapshot end time plus scalar values selected
from the existing result objects. It must not import audio-device or UI code.

- [ ] **Step 2: Document the public contract**

Add a Rolling Analysis section to `docs/python-api.md` covering constructor
arguments, enum values, snapshot fields, window-relative nested timestamps,
reset, mono input, transactional behavior, synchronous/thread ownership,
performance guidance, and no incomplete-window flush. Correct stale version
examples encountered in the directly touched native-loading section.

- [ ] **Step 3: Add concise README discovery**

Add a short paragraph and example link near the Python introduction. Describe
the API as rolling analysis, not native stateful streaming or hard real time.

- [ ] **Step 4: Run complete verification**

Run:

```powershell
$env:PYTHONPATH="$PWD/python"
$env:METASONA_LIBRARY="$PWD/build/py3-none-win_amd64/metasona.dll"
python -m unittest discover -s tests -p test_rolling.py -v
python examples/metasona_rolling.py
python -m compileall -q python examples tests
git diff --check
uv lock --check
```

Then build a wheel using the repository's MSVC developer environment, install
it into a fresh uv Python 3.9 environment, and run
`.github/scripts/smoke_test.py` against that installation.
Expected: every command exits zero and the installed package/native versions
both report 0.2.2.

- [ ] **Step 5: Review changed-line estimate and commit**

Compare the final diff against the spec estimate of 685-990 changed lines and
record any material deviation in the completion report.

```powershell
git add README.md docs/python-api.md examples/metasona_rolling.py
git commit -m "Document rolling analysis workflow"
```
