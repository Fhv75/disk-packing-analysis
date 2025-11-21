"""
Interfaz para análisis y presentación de resultados.

Provee funciones de alto nivel para:
- Imprimir resúmenes en consola
- Mostrar análisis completos con gráficos
- Crear diccionarios de datos para APIs web
"""

from __future__ import annotations
import numpy as np
from .analysis import AnalysisResult


def print_analysis_summary(result: AnalysisResult) -> None:
    """
    Imprime un resumen del análisis en consola.
    
    Args:
        result: Resultado completo del análisis
    """
    config = result.config
    print(f"\n{'='*60}")
    print(f"Análisis de configuración: {config.name or 'Sin nombre'}")
    print(f"{'='*60}")
    print(f"Número de discos: {config.n}")
    print(f"Número de contactos: {len(config.edges)}")
    print(f"Dimensión del rolling space: {result.R.shape[1]}")
    print(f"Perímetro (centros): {result.perimeter_centers:.4f}")
    print(f"Perímetro (discos): {result.perimeter_disks:.4f}")
    print(f"\nAutovalores del Hessiano intrínseco:")
    for i, val in enumerate(result.eigenvalues):
        print(f"  λ[{i}] = {val:.6e}")
    print(f"{'='*60}\n")


def show_complete_analysis(result: AnalysisResult) -> None:
    """
    Muestra análisis completo con todas las gráficas de Matplotlib.
    
    Args:
        result: Resultado completo del análisis
    """
    from .plotting import plot_disks, plot_contact_graph, plot_spectrum
    
    print_analysis_summary(result)
    plot_disks(result.config, show_hull=True)
    plot_contact_graph(result.config, show_normals=True)
    plot_spectrum(result.eigenvalues, result.config.name or "")


def create_dashboard(result: AnalysisResult) -> dict:
    """
    Crea un diccionario con toda la información necesaria para la interfaz web.
    
    Args:
        result: Resultado completo del análisis
        
    Returns:
        Diccionario con datos serializables para JSON
    """
    config = result.config
    
    return {
        'config_name': config.name or 'Sin nombre',
        'n_disks': config.n,
        'n_contacts': len(config.edges),
        'coords': config.coords.tolist(),
        'edges': config.edges,
        'hull_vertices': config.hull_vertices,
        'perimeter_centers': float(result.perimeter_centers),
        'perimeter_disks': float(result.perimeter_disks),
        'contact_matrix': result.A.tolist(),
        'rolling_space_dim': result.R.shape[1],
        'eigenvalues': result.eigenvalues.tolist(),
    }
