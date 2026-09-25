# SPDX-License-Identifier: Apache-2.0
# MetaSona author: Jiahua Zhang, September 2026.

"""Lazy, ABI-checked loading and declaration of the native C interface."""

from __future__ import annotations

import atexit
import ctypes
import ctypes.util
import os
import sys
import threading
from contextlib import ExitStack, suppress
from importlib import resources
from pathlib import Path
from types import TracebackType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from .exceptions import NativeLibraryError

EXPECTED_ABI_VERSION = 1
EXPECTED_NATIVE_VERSION = "0.2.1"
_ENVIRONMENT_VARIABLE = "METASONA_LIBRARY"
_SYSTEM_FALLBACK_ENVIRONMENT_VARIABLE = "METASONA_ALLOW_SYSTEM_LIBRARY"

DoublePointer = ctypes.POINTER(ctypes.c_double)


class NativeLibrary:
    """Configured native library plus safe accessors for static native data."""

    def __init__(self, library: ctypes.CDLL, origin: str) -> None:
        self.lib = library
        self.origin = origin

    @property
    def version(self) -> str:
        raw = self.lib.ms_version_string()
        return raw.decode("utf-8", errors="replace") if raw else "unknown"

    def status_message(self, status: int) -> str:
        raw = self.lib.ms_status_string(status)
        return raw.decode("utf-8", errors="replace") if raw else "unknown native status"

    def bark_axis(self) -> NDArray[np.float64]:
        pointer = self.lib.ms_bark_axis()
        if not pointer:
            raise NativeLibraryError("the native library returned a null Bark-axis pointer")
        return np.ctypeslib.as_array(pointer, shape=(240,)).astype(np.float64, copy=True)

    def third_octave_centres_hz(self) -> NDArray[np.float64]:
        pointer = self.lib.ms_third_octave_centres_hz()
        if not pointer:
            raise NativeLibraryError(
                "the native library returned a null third-octave-frequency pointer"
            )
        return np.ctypeslib.as_array(pointer, shape=(28,)).astype(np.float64, copy=True)


_LOAD_LOCK = threading.Lock()
_NATIVE_LIBRARY: NativeLibrary | None = None
_RESOURCE_CONTEXTS = ExitStack()
_DLL_DIRECTORY_HANDLES: list[Any] = []


def _cleanup_resources() -> None:
    while _DLL_DIRECTORY_HANDLES:
        handle = _DLL_DIRECTORY_HANDLES.pop()
        with suppress(AttributeError, OSError):
            handle.close()
    _RESOURCE_CONTEXTS.close()


atexit.register(_cleanup_resources)


def _configure_signatures(library: ctypes.CDLL) -> None:
    """Declare every public symbol exactly as specified by ABI version 1."""
    library.ms_abi_version.argtypes = []
    library.ms_abi_version.restype = ctypes.c_uint32
    library.ms_version_string.argtypes = []
    library.ms_version_string.restype = ctypes.c_char_p
    library.ms_status_string.argtypes = [ctypes.c_int32]
    library.ms_status_string.restype = ctypes.c_char_p
    library.ms_third_octave_centres_hz.argtypes = []
    library.ms_third_octave_centres_hz.restype = DoublePointer
    library.ms_bark_axis.argtypes = []
    library.ms_bark_axis.restype = DoublePointer

    library.ms_loudness_from_levels.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        ctypes.c_uint32,
        DoublePointer,
        DoublePointer,
        ctypes.c_size_t,
    ]
    library.ms_loudness_from_levels.restype = ctypes.c_int32
    library.ms_loudness_stationary.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.c_uint32,
        DoublePointer,
        DoublePointer,
        ctypes.c_size_t,
    ]
    library.ms_loudness_stationary.restype = ctypes.c_int32
    library.ms_loudness_time_frame_count.argtypes = [
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_loudness_time_frame_count.restype = ctypes.c_int32
    library.ms_loudness_time.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.c_uint32,
        DoublePointer,
        ctypes.c_size_t,
        DoublePointer,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_loudness_time.restype = ctypes.c_int32
    library.ms_roughness_frame_count.argtypes = [
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_roughness_frame_count.restype = ctypes.c_int32
    library.ms_roughness_dw.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        ctypes.c_uint32,
        DoublePointer,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_roughness_dw.restype = ctypes.c_int32
    library.ms_tonality_frame_count.argtypes = [
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_tonality_frame_count.restype = ctypes.c_int32
    library.ms_tonality_aures.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        ctypes.c_uint32,
        ctypes.c_uint32,
        DoublePointer,
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_size_t),
    ]
    library.ms_tonality_aures.restype = ctypes.c_int32
    library.ms_sharpness_din.argtypes = [
        DoublePointer,
        ctypes.c_size_t,
        DoublePointer,
    ]
    library.ms_sharpness_din.restype = ctypes.c_int32

    for name in ("ms_ecma_tonal_frame_count", "ms_ecma_roughness_frame_count"):
        function = getattr(library, name)
        function.argtypes = [ctypes.c_size_t, ctypes.c_uint32, ctypes.POINTER(ctypes.c_size_t)]
        function.restype = ctypes.c_int32
    for name in ("ms_ecma_tonal_analysis", "ms_roughness_ecma"):
        function = getattr(library, name)
        function.argtypes = [DoublePointer, ctypes.c_size_t, ctypes.c_uint32,
                             ctypes.c_uint32, DoublePointer, ctypes.c_size_t,
                             ctypes.POINTER(ctypes.c_size_t)]
        function.restype = ctypes.c_int32


def _library_filenames() -> tuple[str, ...]:
    if sys.platform == "win32":
        return ("metasona.dll", "libmetasona.dll")
    if sys.platform == "darwin":
        return ("libmetasona.dylib", "metasona.dylib")
    return ("libmetasona.so", "libmetasona.so.1", "libmetasona.so.0.2.1")


def _package_candidates() -> list[Any]:
    try:
        native_directory = resources.files("metasona").joinpath("_native")
        children = list(native_directory.iterdir())
    except (FileNotFoundError, ModuleNotFoundError, NotADirectoryError):
        return []
    names = _library_filenames()
    by_name = {child.name: child for child in children if child.is_file()}
    return [by_name[name] for name in names if name in by_name]


def _load_path(path: Path) -> ctypes.CDLL:
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        _DLL_DIRECTORY_HANDLES.append(os.add_dll_directory(str(path.parent)))
    return ctypes.CDLL(str(path))


def _validate_loaded_library(library: ctypes.CDLL, origin: str) -> NativeLibrary:
    """Configure and verify both binary and semantic package compatibility."""
    _configure_signatures(library)
    actual_abi = int(library.ms_abi_version())
    if actual_abi != EXPECTED_ABI_VERSION:
        raise NativeLibraryError(
            f"native ABI mismatch at {origin}: expected {EXPECTED_ABI_VERSION}, found {actual_abi}"
        )
    raw_version = library.ms_version_string()
    actual_version = raw_version.decode("utf-8", errors="replace") if raw_version else "unknown"
    if actual_version != EXPECTED_NATIVE_VERSION:
        raise NativeLibraryError(
            f"native version mismatch at {origin}: expected {EXPECTED_NATIVE_VERSION}, "
            f"found {actual_version}"
        )
    return NativeLibrary(library, origin)


def _attempt_load() -> NativeLibrary:
    attempts: list[str] = []
    explicit = os.environ.get(_ENVIRONMENT_VARIABLE)
    if explicit:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise NativeLibraryError(
                f"{_ENVIRONMENT_VARIABLE} does not name a native library file: {path}"
            )
        candidates: list[tuple[str, Any]] = [(str(path), path)]
    else:
        candidates = [(f"bundled resource {item.name}", item) for item in _package_candidates()]

    for label, candidate in candidates:
        try:
            if isinstance(candidate, Path):
                path = candidate
            else:
                path = _RESOURCE_CONTEXTS.enter_context(resources.as_file(candidate))
            library = _load_path(path)
            return _validate_loaded_library(library, str(path))
        except (AttributeError, OSError, NativeLibraryError) as exc:
            attempts.append(f"{label}: {exc}")

    allow_system = os.environ.get(_SYSTEM_FALLBACK_ENVIRONMENT_VARIABLE) == "1"
    if not explicit and allow_system:
        located = ctypes.util.find_library("metasona")
        if located:
            try:
                library = ctypes.CDLL(located)
                return _validate_loaded_library(library, located)
            except (AttributeError, OSError, NativeLibraryError) as exc:
                attempts.append(f"system library {located}: {exc}")

    if attempts:
        detail = "; ".join(attempts)
    elif allow_system:
        detail = "no bundled or system library was found"
    else:
        detail = (
            "no bundled library was found and system-library lookup is disabled; "
            f"set {_SYSTEM_FALLBACK_ENVIRONMENT_VARIABLE}=1 to opt in"
        )
    raise NativeLibraryError(
        "unable to load metasona's native library. Reinstall a platform wheel or set "
        f"{_ENVIRONMENT_VARIABLE} to a compatible library. Attempts: {detail}"
    )


def get_native_library() -> NativeLibrary:
    """Load and ABI-check the native library on first numerical API use."""
    global _NATIVE_LIBRARY
    if _NATIVE_LIBRARY is None:
        with _LOAD_LOCK:
            if _NATIVE_LIBRARY is None:
                _NATIVE_LIBRARY = _attempt_load()
    return _NATIVE_LIBRARY


def _reset_for_tests(
    exc_type: type[BaseException] | None = None,
    exc: BaseException | None = None,
    traceback: TracebackType | None = None,
) -> None:
    """Clear only the cached library; private and intentionally test-only."""
    del exc_type, exc, traceback
    global _NATIVE_LIBRARY
    _NATIVE_LIBRARY = None
