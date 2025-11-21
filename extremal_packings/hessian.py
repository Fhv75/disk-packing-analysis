from __future__ import annotations
import numpy as np
from .configurations import Configuration
from .constraints import contact_normals


def local_edge_hessian(u: np.ndarray, dist: float) -> np.ndarray:
    """
    H_ij local para una arista (i,j):
    H_ij = (1/dist) * (I - u u^T)
    donde u es normal unitario y dist = ||c_j - c_i||.
    """
    u = u.reshape(2, 1)
    I = np.eye(2)
    P = I - u @ u.T
    return (1.0 / dist) * P


def build_unconstrained_hessian(config: Configuration) -> np.ndarray:
    """
    Ensambla el Hessiano global K(c) en R^{2n}.
    K_ii = sum_j H_ij
    K_ij = -H_ij si (i,j) es contacto
    """
    n = config.n
    coords = config.coords
    K = np.zeros((2 * n, 2 * n), dtype=float)

    normals = contact_normals(config)
    for (i, j), u in normals.items():
        diff = coords[j] - coords[i]
        dist = np.linalg.norm(diff)
        if dist == 0:
            raise ValueError(f"Centros {i} y {j} coinciden (distancia 0)")
        Hij = local_edge_hessian(u, dist)

        # bloques en K
        ii = slice(2 * i, 2 * i + 2)
        jj = slice(2 * j, 2 * j + 2)

        K[ii, ii] += Hij
        K[jj, jj] += Hij
        K[ii, jj] -= Hij
        K[jj, ii] -= Hij

    return K


def project_to_roll(K: np.ndarray, R: np.ndarray) -> np.ndarray:
    """
    Proyección del Hessiano global K a coordenadas del rolling space:
    H = R^T K R
    """
    return R.T @ K @ R


def intrinsic_spectrum(H: np.ndarray) -> np.ndarray:
    """
    Autovalores (ordenados) del Hessiano intrínseco H.
    """
    # Aproximar autovalores pequeños a cero
    eigenvalues = np.linalg.eigvalsh(H)
    eigenvalues[np.abs(eigenvalues) < 1e-12] = 0
    return eigenvalues
    
