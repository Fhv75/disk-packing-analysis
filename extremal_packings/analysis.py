"""
Pipeline de análisis de configuraciones de discos.

Este módulo implementa el flujo completo de análisis, desde la validación
del grafo hasta el cálculo del espectro del Hessiano intrínseco.

Classes:
    AnalysisResult: Contenedor de todos los resultados del análisis.

Functions:
    analyze_configuration: Pipeline completo de análisis.

Pipeline de Análisis
--------------------
1. **Validación**: Verificar conectividad y grados máximos del grafo.
2. **Matriz de Contacto**: Construir A(c) de dimensión m×2n.
3. **Rolling Space**: Calcular base ortonormal R de ker(A).
4. **Geometría**: Perímetros de envolventes convexas.
5. **Hessiano Global**: Ensamblar K(c) en ℝ²ⁿ.
6. **Hessiano Intrínseco**: Proyectar H = R^T K R.
7. **Espectro**: Calcular autovalores ordenados de H.

Complejidad Computacional
--------------------------
Para n discos y m contactos:
- Validación del grafo: O(n + m)
- Matriz de contacto: O(m·n)
- Rolling space (SVD): O(min(m, 2n)²·max(m, 2n))
- Convex hull: O(n log n)
- Hessiano global: O(m·n)
- Hessiano intrínseco: O(d²·n) donde d = dim(Roll)
- Espectro: O(d³)

Total: Dominado por SVD, típicamente O(n³) para grafos densos.

Examples
--------
>>> from extremal_packings import load_configuration, analyze_configuration
>>> config = load_configuration("D5-7")
>>> result = analyze_configuration(config)
>>> print(result.eigenvalues)
[0.0000e+00 6.1803e-01 1.6180e+00]
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .configurations import Configuration
from .contact_graphs import check_graph_validity
from .constraints import (
    build_contact_matrix,
    rolling_space_basis,
)
from .perimeter import (
    perimeter_centers,
    perimeter_disks,
)
from .hessian import (
    build_unconstrained_hessian,
    project_to_roll,
    intrinsic_spectrum,
)


@dataclass
class AnalysisResult:
    """
    Contenedor con todos los resultados del análisis.
    
    Attributes:
        config (Configuration): Configuración analizada.
        A (np.ndarray): Matriz de contacto de dimensión (m, 2n).
        R (np.ndarray): Base ortonormal del rolling space (2n, d).
        K (np.ndarray): Hessiano global no restringido (2n, 2n).
        H (np.ndarray): Hessiano intrínseco proyectado (d, d).
        eigenvalues (np.ndarray): Autovalores ordenados de H.
        perimeter_centers (float): Perímetro de convex hull de centros.
        perimeter_disks (float): Perímetro del cluster de discos.
    
    Properties:
        rolling_dim (int): Dimensión del rolling space.
        is_rigid (bool): True si rolling_dim == 0.
        is_flexible (bool): True si rolling_dim > 0.
        has_negative_eigenvalue (bool): True si min(eigenvalues) < 0.
    
    Examples:
        >>> result = analyze_configuration(config)
        >>> result.rolling_dim
        3
        >>> result.is_flexible
        True
        >>> result.eigenvalues
        array([0.0, 0.618, 1.618])
    """
    config: Configuration

    # Contact matrix y rolling space
    A: np.ndarray
    R: np.ndarray

    # Hessianos
    K: np.ndarray      # Hessiano global no restringido
    H: np.ndarray      # Hessiano proyectado al rolling space

    # Espectro del Hessiano intrínseco
    eigenvalues: np.ndarray

    # Perímetros (centros y discos)
    perimeter_centers: float
    perimeter_disks: float

    @property
    def rolling_dim(self) -> int:
        """
        Dimensión del rolling space.
        
        Returns:
            Número de grados de libertad infinitesimales = 2n - rank(A).
        """
        return self.R.shape[1]
    
    @property
    def is_rigid(self) -> bool:
        """
        Verifica si la configuración es infinitesimalmente rígida.
        
        Returns:
            True si no hay deformaciones (excepto movimientos rígidos).
        
        Notes:
            Una configuración es rígida si rolling_dim ≤ 3 (en el plano).
            Los 3 grados corresponden a traslaciones (2) y rotación (1).
        """
        return self.rolling_dim <= 3
    
    @property
    def is_flexible(self) -> bool:
        """
        Verifica si la configuración tiene flexibilidad infinitesimal.
        
        Returns:
            True si rolling_dim > 3.
        """
        return not self.is_rigid
    
    @property
    def has_negative_eigenvalue(self) -> bool:
        """
        Verifica si el Hessiano intrínseco tiene autovalores negativos.
        
        Returns:
            True si existe λ < -10⁻¹⁰ (considerando tolerancia numérica).
        
        Notes:
            Un autovalor negativo indica que la configuración no es
            un mínimo local del perímetro en el rolling space.
        """
        if len(self.eigenvalues) == 0:
            return False
        return self.eigenvalues[0] < -1e-10


def analyze_configuration(config: Configuration) -> AnalysisResult:
    """
    Pipeline completo de análisis de una configuración.
    
    Ejecuta todas las etapas del análisis en orden:
    1. Validación del grafo de contacto
    2. Construcción de la matriz de contacto A(c)
    3. Cálculo del rolling space (ker A)
    4. Cálculo de perímetros
    5. Construcción del Hessiano global K(c)
    6. Proyección al Hessiano intrínseco H
    7. Análisis espectral de H
    
    Args:
        config: Configuración a analizar.
    
    Returns:
        Objeto AnalysisResult con todos los resultados.
    
    Raises:
        ValueError: Si el grafo no es conexo o tiene vértices con grado > 6.
        ValueError: Si los centros de algún contacto coinciden.
    
    Examples:
        >>> from extremal_packings import Configuration, analyze_configuration
        >>> import numpy as np
        >>> 
        >>> # Cadena de 3 discos
        >>> coords = np.array([[0,0], [2,0], [4,0]])
        >>> edges = [(0,1), (1,2)]
        >>> config = Configuration(coords=coords, edges=edges)
        >>> result = analyze_configuration(config)
        >>> 
        >>> # Inspeccionar resultados
        >>> print(f"Rolling space dim: {result.R.shape[1]}")
        >>> print(f"Eigenvalues: {result.eigenvalues}")
        >>> print(f"Perimeter: {result.perimeter_disks:.4f}")
    
    Notes:
        - Para configuraciones sin contactos (m=0), Roll(c) = ℝ²ⁿ.
        - Autovalores pequeños (|λ| < 10⁻¹⁰) se consideran cero numérico.
        - El orden de los autovalores es no-decreciente.
    
    See Also:
        build_contact_matrix: Construcción de A(c).
        rolling_space_basis: Cálculo de ker(A).
        build_unconstrained_hessian: Construcción de K(c).
        intrinsic_spectrum: Cálculo de autovalores de H.
    """

    # -------------------------------
    # 1. Validación del grafo
    # -------------------------------
    check_graph_validity(config.n, config.edges)

    # -------------------------------
    # 2. Matriz de contacto A(c)
    # -------------------------------
    A = build_contact_matrix(config)

    # -------------------------------
    # 3. Rolling space R = ker(A)
    # -------------------------------
    if len(config.edges) == 0:
        # Caso sin contactos: ker(A) = R^{2n}
        R = np.eye(2 * config.n)
    else:
        R = rolling_space_basis(A)

    # -------------------------------
    # 4. Perímetros
    # -------------------------------
    p_centers = perimeter_centers(config)
    p_disks = perimeter_disks(config)

    # -------------------------------
    # 5. Hessiano global
    # -------------------------------
    K = build_unconstrained_hessian(config)

    # -------------------------------
    # 6. Hessiano proyectado
    # -------------------------------
    H = project_to_roll(K, R)

    # -------------------------------
    # 7. Autovalores
    # -------------------------------
    eigenvalues = intrinsic_spectrum(H)

    # -------------------------------
    # Empaquetar resultados
    # -------------------------------
    return AnalysisResult(
        config=config,
        A=A,
        R=R,
        K=K,
        H=H,
        eigenvalues=eigenvalues,
        perimeter_centers=p_centers,
        perimeter_disks=p_disks,
    )
