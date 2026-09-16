# Pinned PocketFFT C backend

Upstream: https://github.com/mreineck/pocketfft

Commit: `81d171a6d5562e3aaa2c73489b70f564c633ff81` (C/master branch).
The current default upstream branch is C++; this project deliberately uses
the C backend to retain its C11 toolchain.

Files downloaded verbatim from the pinned commit on 2026-09-16:

| File | SHA-256 (upstream LF bytes) |
|---|---|
| pocketfft.c | `4f9f4095080fdba96d7112bcfe6ed66cb5c094077e92f7544602756eccd0f97c` |
| pocketfft.h | `8c12a7b2467a86c3f8818a5ab85cec2a8e90dfb2ff273cb8cf3ad937ea286ba3` |
| LICENSE.md | `c254b2b5fd30ae3c89d36ba061730563b9760fca04f1b746a1b2c3157da2c789` |

No source modifications. CMake prefixes symbols with `ms_pocketfft_`, sets
hidden visibility, and builds this third-party translation unit separately
from MetaSona's warnings-as-errors policy. The reused metric kernels call
PocketFFT directly and reuse their per-call plans and buffers across
frames/channels. Plans are never shared between calls.

License: BSD-3-Clause, retained in `LICENSE.md` and copied to
`../../../LICENSES/BSD-3-Clause-PocketFFT.txt` (distribution license directory).
