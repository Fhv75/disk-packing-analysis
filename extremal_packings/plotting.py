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


# Paleta de colores consistente
COLORS = {
    'disk_edge': '#2E3A59',      # Azul oscuro para bordes de discos
    'disk_fill': '#E8EAF6',      # Azul muy claro para relleno
    'center': '#1A1A1A',         # Negro para centros
    'contact_tangent': '#9C27B0', # Púrpura para tangentes
    'hull': '#2E3A59',           # Azul oscuro para envolvente
    'label': '#1A1A1A',          # Negro para etiquetas
    'spectrum_bar': '#64B5F6',   # Azul para barras de espectro
}

# Colores pastel para las aristas (cada arista tendrá un color diferente)
EDGE_COLORS = [
    '#77B6EA',  # Azul medio
    '#FFB577',  # Naranja melocotón
    '#88D8B0',  # Verde menta medio
    '#FF9AA2',  # Rosa coral
    '#C38EC7',  # Lavanda medio
    '#D4A5A5',  # Café claro
    '#FF9EC5',  # Rosa medio
    '#B0B0B0',  # Gris medio
    '#FFD97D',  # Amarillo dorado
]


def plot_disks(config: Configuration, ax=None, radius: float = 1.0, 
               show_labels: bool = True, show_centers: bool = True, 
               show_hull: bool = True, show_contacts: bool = True,
               show_tangents: bool = False) -> None:
    """
    Visualiza los discos completos con sus centros y opcionalmente la envolvente convexa.
    
    Args:
        config: Configuración de discos a visualizar
        ax: Axes de matplotlib (si es None, crea una nueva figura)
        radius: Radio de los discos (default: 1.0)
        show_labels: Si se muestran las etiquetas de los discos
        show_centers: Si se muestran los centros marcados
        show_hull: Si se muestra la envolvente convexa DE LOS DISCOS
        show_contacts: Si se muestran las aristas de contacto
        show_tangents: Si se muestran las tangentes de contacto
    """
    coords = config.coords
    
    # Crear figura solo si no se proporciona ax
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
        standalone = True
    else:
        standalone = False
    
    # Establecer fondo blanco
    if standalone:
        fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Dibujar discos PRIMERO
    for i, (xi, yi) in enumerate(coords):
        circle = Circle((xi, yi), radius, 
                       fill=True, 
                       facecolor=COLORS['disk_fill'],
                       edgecolor=COLORS['disk_edge'], 
                       linewidth=2,
                       alpha=0.7,
                       zorder=2)
        ax.add_patch(circle)
    
    # Dibujar aristas de contacto
    if show_contacts:
        for idx, (i, j) in enumerate(config.edges):
            edge_color = EDGE_COLORS[idx % len(EDGE_COLORS)]
            xs = [coords[i, 0], coords[j, 0]]
            ys = [coords[i, 1], coords[j, 1]]
            ax.plot(xs, ys, color=edge_color, 
                   linewidth=2.5, zorder=3, alpha=0.9)
    
    # Dibujar tangentes de contacto
    if show_tangents:
        for idx, (i, j) in enumerate(config.edges):
            ci, cj = coords[i], coords[j]
            diff = cj - ci
            dist = np.linalg.norm(diff)
            if dist > 1e-10:
                contact_point = ci + (diff / dist) * radius
                tangent = np.array([-diff[1], diff[0]]) / dist
                t_len = 0.5
                p1 = contact_point - tangent * t_len
                p2 = contact_point + tangent * t_len
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], 
                       color=COLORS['contact_tangent'], 
                       linewidth=2, zorder=4, alpha=0.9)
    
    # Dibujar centros
    if show_centers:
        ax.scatter(coords[:, 0], coords[:, 1], 
                  color=COLORS['center'], s=50, zorder=5,
                  edgecolors='white', linewidths=1)
    
    # Dibujar etiquetas FUERA de los centros
    if show_labels:
        for i, (xi, yi) in enumerate(coords):
            centroid = coords.mean(axis=0)
            direction = np.array([xi, yi]) - centroid
            direction_norm = np.linalg.norm(direction)
            
            if direction_norm > 1e-10:
                offset_dir = direction / direction_norm * 0.2
                text_x = xi + offset_dir[0]
                text_y = yi + offset_dir[1]
            else:
                text_x = xi
                text_y = yi + 0.2
            
            ax.text(text_x, text_y, str(i), ha="center", va="center", 
                   fontsize=11, fontweight='bold', 
                   color=COLORS['label'], zorder=6,
                   bbox=dict(boxstyle='round,pad=0.3', 
                           facecolor='white', 
                           edgecolor='none',
                           alpha=0.8))
    
    # Dibujar envolvente convexa DE LOS DISCOS (no solo centros)
    if show_hull:
        from .perimeter import compute_disk_hull_geometry
        
        try:
            # Obtener la geometría completa del hull de discos
            hull_geom = compute_disk_hull_geometry(config, radius)
            
            if hull_geom:
                # Dibujar segmentos tangentes
                if 'tangent_segments' in hull_geom:
                    for segment in hull_geom['tangent_segments']:
                        start = segment['start']
                        end = segment['end']
                        ax.plot([start[0], end[0]], [start[1], end[1]], 
                               color=COLORS['hull'], linewidth=4, 
                               linestyle='-', zorder=10, alpha=1.0)
                
                # Dibujar arcos circulares
                if 'arcs' in hull_geom:
                    for arc in hull_geom['arcs']:
                        arc_x, arc_y = _create_arc_points(
                            arc['center'], 
                            arc['radius'], 
                            arc['angle_start'], 
                            arc['angle_end']
                        )
                        ax.plot(arc_x, arc_y, 
                               color=COLORS['hull'], linewidth=4, 
                               linestyle='-', zorder=10, alpha=1.0)
                
                print("DEBUG: Hull de discos dibujado exitosamente")
            else:
                print("DEBUG: hull_geom es None")
                
        except Exception as e:
            print(f"ERROR al dibujar hull de discos: {e}")
            import traceback
            traceback.print_exc()
    
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("y", fontsize=12)
    ax.set_title(f"{config.name} (Discos)", fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    if standalone:
        plt.tight_layout()
        plt.show()


def _create_arc_points(center, radius, angle_start, angle_end, num_points=100):
    """
    Crea puntos para dibujar un arco circular EXTERIOR del hull.
    Siempre toma el arco < 180° (el más corto).
    
    Args:
        center: Centro del círculo [x, y] o (x, y)
        radius: Radio del círculo
        angle_start: Ángulo inicial en grados
        angle_end: Ángulo final en grados
        num_points: Número de puntos del arco
    
    Returns:
        (x_points, y_points): Arrays de coordenadas del arco
    """
    # Convertir a array si es necesario
    if isinstance(center, (list, tuple)):
        center = np.array(center)
    
    # Convertir ángulos a radianes
    angle_start_rad = np.radians(angle_start)
    angle_end_rad = np.radians(angle_end)
    
    # Calcular diferencia de ángulo
    angle_diff = angle_end_rad - angle_start_rad
    
    # Si es negativo, agregar 2π para ir en sentido positivo
    if angle_diff < 0:
        angle_diff += 2 * np.pi
    
    # CLAVE: Si el arco es > 180°, invertir para tomar el arco complementario
    # (el más corto, que es el EXTERIOR del cluster)
    if angle_diff > np.pi:
        # Invertir: ir desde angle_end hasta angle_start + 2π
        final_start_rad = angle_end_rad
        final_end_rad = angle_start_rad + 2 * np.pi
        final_diff = 2 * np.pi - angle_diff
    else:
        # El arco es < 180°, usarlo directamente
        final_start_rad = angle_start_rad
        final_end_rad = angle_end_rad
        final_diff = angle_diff
    
    # Generar ángulos uniformemente espaciados
    angles = np.linspace(final_start_rad, final_start_rad + final_diff, num_points)
    
    # Calcular puntos del arco
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    
    return x, y

def plot_contact_graph(config: Configuration, ax=None, show_normals: bool = False) -> None:
    """
    Visualiza el grafo de contacto con opciones para mostrar vectores normales.
    
    Args:
        config: Configuración de discos
        ax: Axes de matplotlib (si es None, crea una nueva figura)
        show_normals: Si se muestran los vectores normales de contacto
    """
    coords = config.coords
    
    # Crear figura solo si no se proporciona ax
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 10))
        standalone = True
    else:
        standalone = False
    
    # Establecer fondo blanco
    if standalone:
        fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Dibujar aristas de contacto (cada una con color diferente)
    for idx, (i, j) in enumerate(config.edges):
        edge_color = EDGE_COLORS[idx % len(EDGE_COLORS)]
        xs = [coords[i, 0], coords[j, 0]]
        ys = [coords[i, 1], coords[j, 1]]
        ax.plot(xs, ys, color=edge_color, 
               linewidth=3, zorder=2, alpha=0.9)
        
        # Dibujar vectores normales si se solicita
        if show_normals:
            diff = coords[j] - coords[i]
            dist = np.linalg.norm(diff)
            if dist > 0:
                u = diff / dist
                mid = (coords[i] + coords[j]) / 2
                ax.arrow(mid[0], mid[1], u[0] * 0.4, u[1] * 0.4,
                        head_width=0.15, head_length=0.15, 
                        fc=COLORS['contact_tangent'], 
                        ec=COLORS['contact_tangent'],
                        linewidth=2, zorder=3)
    
    # Dibujar nodos (centros)
    ax.scatter(coords[:, 0], coords[:, 1], 
              s=400, c=COLORS['center'], zorder=4,
              edgecolors='white', linewidths=2)
    
    # Dibujar etiquetas DENTRO de los nodos (en blanco sobre negro)
    for k, (xi, yi) in enumerate(coords):
        ax.text(xi, yi, str(k), ha="center", va="center", 
               fontsize=13, fontweight='bold', color='white', zorder=5)
    
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("y", fontsize=12)
    ax.set_title(f"{config.name} (Centros)", 
                fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.2, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Solo mostrar si es standalone
    if standalone:
        plt.tight_layout()
        plt.show()


def plot_spectrum(eigenvalues: np.ndarray, ax=None, config_name: str = "") -> None:
    """
    Visualiza el espectro (autovalores) del Hessiano intrínseco.
    
    Args:
        eigenvalues: Array de autovalores ordenados
        ax: Axes de matplotlib (si es None, crea una nueva figura)
        config_name: Nombre de la configuración para el título
    """
    # Crear figura solo si no se proporciona ax
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 7))
        standalone = True
    else:
        standalone = False
    
    # Establecer fondo blanco
    if standalone:
        fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    x = np.arange(len(eigenvalues))
    
    # Colores diferentes para valores positivos y negativos
    colors = [COLORS['spectrum_bar'] if val >= 0 else COLORS['contact_tangent'] 
              for val in eigenvalues]
    
    bars = ax.bar(x, eigenvalues, color=colors, edgecolor=COLORS['disk_edge'], 
                  linewidth=1.5, alpha=0.8)
    
    # Línea de referencia en cero
    ax.axhline(y=0, color=COLORS['disk_edge'], linestyle='--', 
              linewidth=2, alpha=0.7)
    
    ax.set_xlabel("Índice del autovalor", fontsize=12, fontweight='bold')
    ax.set_ylabel("Autovalor", fontsize=12, fontweight='bold')
    ax.set_title(f"Espectro del Hessiano Intrínseco - {config_name}", 
                fontsize=14, fontweight='bold', pad=15)
    ax.grid(True, alpha=0.2, axis='y', linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Añadir valores en las barras (solo para valores significativos)
    for i, (val, bar) in enumerate(zip(eigenvalues, bars)):
        if abs(val) > 1e-10:  # Solo mostrar si no es prácticamente cero
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{val:.2e}',
                   ha='center', va='bottom' if val > 0 else 'top',
                   fontsize=9, fontweight='bold')
    
    # Solo mostrar si es standalone
    if standalone:
        plt.tight_layout()
        plt.show()


def plot_complete_analysis(config: Configuration, result=None, figsize=(16, 12)):
    """
    Crea una visualización completa con 4 subplots:
    1. Discos con hull
    2. Grafo de contacto
    3. Discos con tangentes
    4. Espectro
    
    Args:
        config: Configuración a visualizar
        result: Resultado del análisis (si None, se calcula)
        figsize: Tamaño de la figura completa
    """
    from .analysis import analyze_configuration
    
    if result is None:
        result = analyze_configuration(config)
    
    fig = plt.figure(figsize=figsize)
    fig.patch.set_facecolor('white')
    
    # Subplot 1: Discos con hull y contactos
    ax1 = plt.subplot(2, 2, 1)
    plot_disks(config, ax=ax1, show_hull=True, show_tangents=False)
    
    # Subplot 2: Grafo de contacto
    ax2 = plt.subplot(2, 2, 2)
    plot_contact_graph(config, ax=ax2, show_normals=True)
    
    # Subplot 3: Discos con tangentes
    ax3 = plt.subplot(2, 2, 3)
    plot_disks(config, ax=ax3, show_hull=False, show_tangents=True)
    
    # Subplot 4: Espectro
    ax4 = plt.subplot(2, 2, 4)
    plot_spectrum(result.eigenvalues, ax=ax4, config_name=config.name)
    
    plt.suptitle(f"Análisis Completo: {config.name}", 
                fontsize=16, fontweight='bold', y=0.995)
    plt.tight_layout()
    plt.show()