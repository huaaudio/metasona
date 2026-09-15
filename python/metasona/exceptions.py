# SPDX-License-Identifier: Apache-2.0

"""Exception hierarchy for :mod:`metasona`."""

from __future__ import annotations


class MetaSonaError(Exception):
    """Base class for all package-specific exceptions."""


class MetaSonaValidationError(MetaSonaError, ValueError):
    """Raised before native execution when an input violates the public contract."""


class NativeLibraryError(MetaSonaError, OSError):
    """Raised when the bundled native library cannot be loaded or has the wrong ABI."""


class NativeCallError(MetaSonaError, RuntimeError):
    """Raised when a validated native operation returns an error status."""

    def __init__(self, operation: str, status: int, message: str) -> None:
        self.operation = operation
        self.status = status
        self.native_message = message
        super().__init__(f"{operation} failed with native status {status}: {message}")
