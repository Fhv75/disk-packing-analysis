"""
Módulo de visualización con Matplotlib para análisis de configuraciones de discos.

Este módulo provee funciones para visualización estática usando Matplotlib.
Se utiliza principalmente para análisis en notebooks o scripts de línea de comandos.
"""

from __future__ import annotations
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from .configurations import Configuration
from .perimeter import compute_hull


def plot_disks(config: Configuration, radius: float = 1.0, show_labels: bool = True, 
               show_centers: bool = True, show_hull: bool = False) -> None:
    """
    Visualiza los discos completos con sus centros y opcionalmente la envolvente convexa.
    
    Args:
        config: Configuración de discos a visualizar
        radius: Radio de los discos (default: 1.0)
        show_labels: Si se muestran las etiquetas de los discos
        show_centers: Si se muestran los centros marcados
        show_hull: Si se muestra la envolvente convexa
    """
    coords = config.coords
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Dibujar discos
    for i, (xi, yi) in enumerate(coords):
        circle = Circle((xi, yi), radius, fill=False, edgecolor='blue', linewidth=1.5)
        ax.add_patch(circle)
        if show_labels:
            ax.text(xi, yi, str(i), ha="center", va="center", fontsize=10, fontweight='bold')
    
    # Dibujar centros
    if show_centers:
        ax.scatter(coords[:, 0], coords[:, 1], color='red', s=30, zorder=5)
    
    # Dibujar envolvente convexa
    if show_hull:
        hull = config.hull_vertices if config.hull_vertices is not None else compute_hull(config)
        if len(hull) >= 2:
            for k in range(len(hull)):
                i = hull[k]
                j = hull[(k + 1) % len(hull)]
                xs = [coords[i, 0], coords[j, 0]]
                ys = [coords[i, 1], coords[j, 1]]
                ax.plot(xs, ys, color='green', linewidth=2, linestyle='--')
    
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(config.name or "Disk Configuration")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_contact_graph(config: Configuration, show_normals: bool = False) -> None:
    """
    Visualiza el grafo de contacto con opciones para mostrar vectores normales.
    
    Args:
        config: Configuración de discos
        show_normals: Si se muestran los vectores normales de contacto
    """
    coords = config.coords
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(coords[:, 0], coords[:, 1], s=100, c='red', zorder=5)
    
    # Dibujar aristas de contacto
    for (i, j) in config.edges:
        xs = [coords[i, 0], coords[j, 0]]
        ys = [coords[i, 1], coords[j, 1]]
        ax.plot(xs, ys, 'b-', linewidth=2)
        
        # Dibujar vectores normales si se solicita
        if show_normals:
            diff = coords[j] - coords[i]
            dist = np.linalg.norm(diff)
            if dist > 0:
                u = diff / dist
                mid = (coords[i] + coords[j]) / 2
                ax.arrow(mid[0], mid[1], u[0] * 0.3, u[1] * 0.3,
                        head_width=0.1, head_length=0.1, fc='green', ec='green')
    
    for k, (xi, yi) in enumerate(coords):
        ax.text(xi, yi, str(k), ha="center", va="center", fontsize=12, 
                fontweight='bold', color='white')
    
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(config.name or "Contact graph")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def plot_spectrum(eigenvalues: np.ndarray, config_name: str = "") -> None:
    """
    Visualiza el espectro (autovalores) del Hessiano intrínseco.
    
    Args:
        eigenvalues: Array de autovalores ordenados
        config_name: Nombre de la configuración para el título
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(eigenvalues))
    ax.bar(x, eigenvalues, color='steelblue', edgecolor='black')
    ax.axhline(y=0, color='r', linestyle='--', linewidth=1)
    ax.set_xlabel("Índice del autovalor")
    ax.set_ylabel("Autovalor")
    ax.set_title(f"Espectro del Hessiano Intrínseco {config_name}")
    ax.grid(True, alpha=0.3, axis='y')
    
    # Añadir valores en las barras
    for i, val in enumerate(eigenvalues):
        ax.text(i, val, f'{val:.2e}', ha='center', 
                va='bottom' if val > 0 else 'top', fontsize=8)
    
    plt.tight_layout()
    plt.show()