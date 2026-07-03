"""GPU trace-of-Frobenius tables (optional backend for ecfp._trace_table).

Port of the user's Colab kernel (2023): build the evaluation matrix
X[j][x] = x^3 + f_j x + g_j mod p in dense blocks and compute the quadratic
character of every entry SIMULTANEOUSLY by elementwise square-and-multiply,
chi(X) = X^((p-1)/2) mod p -- O(log p) dense passes, which is what GPUs like --
then row-sum into all p traces at once.

Differences from the Colab original, so this is a drop-in for the pipeline:
  * rows are the pipeline's canonical models j_to_fg(j) (f = -3j(j-1728),
    g = 2j(j-1728)^2), not the y^2 = x^3 - 3t x + 2t family, so the SIGNS
    match what trfr_to_models expects;
  * int64 throughout (the float path is only exact while p^2 < 2^24, and MPS
    has no float64); exact for p < 3e9;
  * row blocks, so memory stays ~block_elems * 8 bytes instead of p^2 * 8.

Usage (e.g. at the top of a factory notebook):

    import trace_gpu
    trace_gpu.enable()          # CUDA if present, else MPS, else stays numpy

After enable(), every trfr_to_js / trfr_to_models / class scan anywhere in the
pipeline uses the GPU table for p >= min_p; the numpy path remains the
fallback and the source of truth (results are identical).
"""

import numpy as np


def _pick_device():
    import torch
    if torch.cuda.is_available():
        return 'cuda'
    if torch.backends.mps.is_available():
        return 'mps'
    return None


def trace_table_torch(p: int, device: str = None, block_elems: int = 1 << 24):
    """tr[j] = trace of Frobenius of j_to_fg(j) over F_p, as a numpy int64 array.
    Same contract as ecfp._trace_table (entries at j = 0, 1728 are unused).

    The Colab original computed the character by elementwise powering
    X^((p-1)/2) -- ideal on CUDA float tensors, but ~2 log2(p) dense int64
    passes, and Metal's int64 throughput makes that LOSE to numpy.  This uses
    the same dense-block framing with a precomputed chi lookup table instead:
    3 arithmetic passes + 1 gather + 1 reduction per block."""
    import torch
    if device is None:
        device = _pick_device()
        if device is None:
            raise RuntimeError('no GPU backend (cuda/mps) available')
    dev = torch.device(device)
    x = torch.arange(p, dtype=torch.int64, device=dev)
    chi = torch.full((p,), -1, dtype=torch.int8, device=dev)
    chi[(x * x) % p] = 1
    chi[0] = 0
    cubes = (x * x % p) * x % p
    t = (x - 1728) % p
    f = (-3 * x % p) * t % p                   # j_to_fg, vectorized
    g = (2 * x % p) * (t * t % p) % p
    tr = torch.empty(p, dtype=torch.int64, device=dev)
    B = max(1, block_elems // p)
    for s in range(0, p, B):
        E = min(s + B, p)
        X = (cubes.unsqueeze(0) + f[s:E].unsqueeze(1) * x.unsqueeze(0)
             + g[s:E].unsqueeze(1)) % p
        tr[s:E] = -chi[X].sum(dim=1, dtype=torch.int64)
    return tr.cpu().numpy()


def enable(min_p: int = 16384, device: str = None, block_elems: int = 1 << 24) -> str:
    """Plug the GPU table into ecfp; below min_p the numpy path keeps winning
    (kernel-launch overhead -- measured crossover on Apple MPS is ~16k; CUDA
    machines can likely lower this).  Returns the device used, or 'numpy'."""
    import ecfp
    if device is None:
        device = _pick_device()
    if device is None:
        return 'numpy'
    def impl(p):
        if p < min_p:
            return None                        # fall through to numpy
        return trace_table_torch(p, device, block_elems)
    ecfp.TRACE_TABLE_IMPL = impl
    ecfp._trace_table.cache_clear()
    return device
