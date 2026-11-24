"""Nelder-Mead simplex minimisation.

The only optimiser pydock needs is the final polish of a GA solution, and
pulling in SciPy for one derivative-free minimiser costs more than a hundred
megabytes in a packaged application.  This is the adaptive-parameter variant
of Gao & Han, Comput. Optim. Appl. 2012, 51, 259 -- the same algorithm, and
the same parameter choice, that ``scipy.optimize.minimize(method="Nelder-Mead",
options={"adaptive": True})`` uses, so results match within tolerance.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Result:
    x: np.ndarray
    fun: float
    nit: int
    nfev: int
    success: bool


def nelder_mead(func, x0, maxiter: int | None = None, maxfev: int | None = None,
                xatol: float = 1e-4, fatol: float = 1e-4,
                initial_step: float = 0.05, zero_step: float = 0.00025,
                adaptive: bool = True) -> Result:
    """Minimise `func` from `x0` without derivatives.

    Parameters follow SciPy's so this is a drop-in replacement:
    `xatol`/`fatol` are the simplex size and spread at which to stop.
    """
    x0 = np.asarray(x0, dtype=np.float64).ravel()
    n = x0.size
    if maxiter is None:
        maxiter = n * 200
    if maxfev is None:
        maxfev = n * 200

    if adaptive and n > 0:
        # Parameters scaled with dimension: markedly more robust above ~10
        # variables, which every docking chromosome is.
        alpha, gamma, rho, sigma = 1.0, 1.0 + 2.0 / n, 0.75 - 1.0 / (2 * n), 1.0 - 1.0 / n
    else:
        alpha, gamma, rho, sigma = 1.0, 2.0, 0.5, 0.5

    # Initial simplex: SciPy's rule, a 5% step per coordinate.
    sim = np.empty((n + 1, n), dtype=np.float64)
    sim[0] = x0
    for k in range(n):
        y = x0.copy()
        y[k] = y[k] * (1 + initial_step) if y[k] != 0 else zero_step
        sim[k + 1] = y

    nfev = 0

    def f(x):
        nonlocal nfev
        nfev += 1
        v = float(func(x))
        return v if np.isfinite(v) else np.inf

    fsim = np.array([f(s) for s in sim])
    order = np.argsort(fsim)
    sim, fsim = sim[order], fsim[order]

    it = 0
    while it < maxiter and nfev < maxfev:
        # Converged when the simplex is both small and flat.
        if (np.max(np.abs(sim[1:] - sim[0])) <= xatol
                and np.max(np.abs(fsim[1:] - fsim[0])) <= fatol):
            break

        centroid = sim[:-1].mean(axis=0)
        worst = sim[-1]

        xr = centroid + alpha * (centroid - worst)
        fr = f(xr)

        if fr < fsim[0]:
            xe = centroid + gamma * (xr - centroid)
            fe = f(xe)
            sim[-1], fsim[-1] = (xe, fe) if fe < fr else (xr, fr)
        elif fr < fsim[-2]:
            sim[-1], fsim[-1] = xr, fr
        else:
            if fr < fsim[-1]:
                xc = centroid + rho * (xr - centroid)       # outside contraction
                fc = f(xc)
                shrink = fc > fr
            else:
                xc = centroid - rho * (centroid - worst)    # inside contraction
                fc = f(xc)
                shrink = fc >= fsim[-1]
            if shrink:
                sim[1:] = sim[0] + sigma * (sim[1:] - sim[0])
                for k in range(1, n + 1):
                    fsim[k] = f(sim[k])
            else:
                sim[-1], fsim[-1] = xc, fc

        order = np.argsort(fsim)
        sim, fsim = sim[order], fsim[order]
        it += 1

    return Result(x=sim[0].copy(), fun=float(fsim[0]), nit=it, nfev=nfev,
                  success=it < maxiter and nfev <= maxfev)
