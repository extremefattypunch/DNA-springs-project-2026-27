"""Minimal dependency-free MRC / MRCS reader and writer.

Written instead of installing `mrcfile` because the only Python interpreter on this
cluster with a scientific stack is CryoSPARC's vendored worker env, which must not be
modified, and /n/holylabs is at 97.5% of its quota so a new conda env is not an option.

Volumes and stacks are returned as float32 arrays indexed [z, y, x] (or [n, y, x]).
"""

import numpy as np
import struct

_MODE_DTYPE = {0: np.int8, 1: np.int16, 2: np.float32, 6: np.uint16, 12: np.float16}


def read_mrc(path):
    """Return (data[z,y,x] float32, apix) for an MRC/MRCS file."""
    with open(path, "rb") as f:
        h = f.read(1024)
        nx, ny, nz, mode = struct.unpack("<4i", h[0:16])
        mx, my, mz = struct.unpack("<3i", h[28:40])
        xlen, ylen, zlen = struct.unpack("<3f", h[40:52])
        mapc, mapr, maps_ = struct.unpack("<3i", h[64:76])
        nsymbt = struct.unpack("<i", h[92:96])[0]
        if mode not in _MODE_DTYPE:
            raise ValueError(f"{path}: unsupported MRC mode {mode}")
        f.seek(1024 + nsymbt)
        data = np.fromfile(f, dtype=_MODE_DTYPE[mode], count=nx * ny * nz)
    if data.size != nx * ny * nz:
        raise ValueError(f"{path}: truncated ({data.size} of {nx*ny*nz} voxels)")
    data = data.reshape(nz, ny, nx).astype(np.float32)
    if (mapc, mapr, maps_) != (1, 2, 3):
        raise ValueError(f"{path}: unsupported axis order {(mapc, mapr, maps_)}")
    apix = xlen / mx if mx else 1.0
    return data, float(apix)


def write_mrc(path, data, apix, is_stack=False):
    """Write a float32 MRC volume (or MRCS stack when is_stack=True)."""
    data = np.ascontiguousarray(data, dtype=np.float32)
    if data.ndim != 3:
        raise ValueError("data must be 3D [z,y,x] or [n,y,x]")
    nz, ny, nx = data.shape
    h = bytearray(1024)
    struct.pack_into("<4i", h, 0, nx, ny, nz, 2)          # dims + mode 2 (float32)
    struct.pack_into("<3i", h, 16, 0, 0, 0)               # nxstart/nystart/nzstart
    # For an image stack the cell along z is one pixel per image, not a sampled axis.
    mz = 1 if is_stack else nz
    struct.pack_into("<3i", h, 28, nx, ny, mz)            # mx, my, mz
    struct.pack_into("<3f", h, 40, nx * apix, ny * apix, mz * apix)
    struct.pack_into("<3f", h, 52, 90.0, 90.0, 90.0)      # cell angles
    struct.pack_into("<3i", h, 64, 1, 2, 3)               # mapc, mapr, maps
    struct.pack_into("<3f", h, 76, float(data.min()), float(data.max()), float(data.mean()))
    struct.pack_into("<2i", h, 88, 0 if is_stack else 1, 0)   # ispg, nsymbt
    struct.pack_into("<3f", h, 196, 0.0, 0.0, 0.0)        # origin
    struct.pack_into("<4s", h, 208, b"MAP ")
    struct.pack_into("<4B", h, 212, 0x44, 0x44, 0x00, 0x00)   # little-endian machine stamp
    struct.pack_into("<f", h, 216, float(data.std()))     # rms
    struct.pack_into("<i", h, 220, 0)                     # nlabl
    with open(path, "wb") as f:
        f.write(bytes(h))
        data.tofile(f)


def write_mrcs(path, stack, apix):
    """Write a particle stack [n, y, x] as an MRCS image stack."""
    write_mrc(path, stack, apix, is_stack=True)
