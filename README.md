# Extremal Packings - Análisis de Discos Tangentes

Paquete Python para análisis geométrico y espectral de configuraciones de discos unitarios tangentes en el plano.

## 🎯 Características

- **Catálogo de configuraciones predefinidas de 3 a 6 discos**
- **Análisis variacional**: Matriz de contacto, rolling space, proyección del gradiente, Hessiano intrínseco
- **Visualización**: Gráficos de discos, grafos de contacto, espectros
- **API intuitiva**: Funciones de alto nivel para análisis rápido

## 📦 Instalación

### Desde PyPI (cuando esté publicado)
```bash
pip install extremal-packings
```

### Desde el repositorio
```bash
git clone https://github.com/fhenr/disk-packing-analysis.git
cd disk-packing-analysis
pip install -e .
```

## 🚀 Uso Rápido

```python
from extremal_packings import load_configuration, analyze_configuration

# Cargar configuración del catálogo
config = load_configuration("D5-7")  # Pentágono regular

# Análisis completo
result = analyze_configuration(config)

# Ver resultados
print(f"Rolling space dimension: {result.R.shape[1]}")
print(f"Eigenvalues: {result.eigenvalues}")
print(f"Perimeter: {result.perimeter_disks:.4f}")
```

## 📖 Conceptos Clave

### Configuración
n discos unitarios con centros 
$c_1, ..., c_n \in \mathbb{R}^2$ y grafo de contacto $G$ donde dos discos se tocan si $||c_j - c_i|| = 2$

### Matriz de Contacto $A(c)$
Matriz $m×2n$ que codifica restricciones de contacto. Cada fila representa un contacto $(i,j)$.

### Rolling Space
$\text{Roll}(c) = \ker(A(c)) ⊆ ℝ²ⁿ$, Espacio de deformaciones infinitesimales que preservan contactos.

### Hessiano Intrínseco
$H = R^T K(c) R$, Proyección del Hessiano del perímetro al rolling space. Sus autovalores indican estabilidad.

## 📖 Documentación

- **[API Reference](docs/api.md)** - Referencia completa de la API
- **[Documentación Detallada](docs/DETAILED_DOCUMENTATION.md)** - Guía exhaustiva con fundamentos matemáticos
- **[Ejemplos Básicos](examples/basic_usage.py)** - Casos de uso comunes
- **[Ejemplos Avanzados](examples/advanced_usage.py)** - Análisis detallado paso a paso

## 📊 Ejemplos

Ver [`examples/basic_usage.py`](examples/basic_usage.py) para casos de uso completos:

1. Cargar del catálogo
2. Crear configuraciones personalizadas
3. Análisis completo
4. Visualización
5. Análisis comparativo

## 🔬 Estructura del Proyecto

```
extremal_packings/
├── __init__.py          # API pública
├── configurations.py    # Clase Configuration
├── analysis.py          # Pipeline de análisis
├── constraints.py       # Matriz de contacto y rolling space
├── hessian.py          # Cálculo del Hessiano
├── perimeter.py        # Perímetros y convex hull
├── catalog.py          # Catálogo de configuraciones
├── plotting.py         # Visualización
└── interface.py        # Funciones de alto nivel
```

## 🧪 Testing

```bash
pytest tests/
```

## 📝 Licencia

MIT License - Ver [LICENSE](LICENSE) para detalles.

## 👤 Autor

**Fabián Andrés Henry Vilaxa**
**Jose Ayala Hoffman**

## 📈 Roadmap

- [ ] Publicar en PyPI
- [ ] Añadir más configuraciones al catálogo
- [ ] Implementar visualización interactiva
- [ ] Soporte para análisis batch
- [ ] Documentación completa con Sphinx