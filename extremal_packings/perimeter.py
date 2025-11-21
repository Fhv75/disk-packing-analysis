from __future__ import annotations
from typing import List
import numpy as np
from .configurations import Configuration


def _cross(o: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Producto cruzado 2D (o->a) x (o->b)."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def compute_hull(config: Configuration) -> List[int]:
    """
    Envolvente convexa de los centros (monotone chain).
    Devuelve índices de los vértices del hull en orden CCW.
    """
    pts = config.coords
    n = pts.shape[0]
    if n <= 1:
        return list(range(n))
    
    if n == 2:
        return [0, 1]

    # ordenar por (x,y)
    idx_sorted = sorted(range(n), key=lambda i: (pts[i, 0], pts[i, 1]))
    
    lower: List[int] = []
    for i in idx_sorted:
        while len(lower) >= 2 and _cross(
            pts[lower[-2]], pts[lower[-1]], pts[i]
        ) <= 0:
            lower.pop()
        lower.append(i)

    upper: List[int] = []
    for i in reversed(idx_sorted):
        while len(upper) >= 2 and _cross(
            pts[upper[-2]], pts[upper[-1]], pts[i]
        ) <= 0:
            upper.pop()
        upper.append(i)

    # Remover el último punto de cada lista porque se duplica
    hull = lower[:-1] + upper[:-1]
    
    # Quitar duplicados preservando orden
    seen = set()
    unique_hull: List[int] = []
    for i in hull:
        if i not in seen:
            seen.add(i)
            unique_hull.append(i)
    
    return unique_hull


def perimeter_centers(config: Configuration) -> float:
    """Perímetro de la envolvente convexa de los centros."""
    coords = config.coords
    hull = config.hull_vertices if config.hull_vertices is not None else compute_hull(config)
    if len(hull) <= 1:
        return 0.0
    per = 0.0
    for k in range(len(hull)):
        i = hull[k]
        j = hull[(k + 1) % len(hull)]
        per += float(np.linalg.norm(coords[j] - coords[i]))
    return per


def perimeter_disks(config: Configuration, radius: float = 1.0) -> float:
    """
    Perímetro del cluster de discos: Per(center hull) + 2πr.
    Por defecto r = 1 (discos unitarios).
    """
    return perimeter_centers(config) + 2.0 * np.pi * radius
