# Extremal Packings - Análisis de Discos Tangentes

Paquete Python para análisis geométrico y espectral de configuraciones de discos unitarios tangentes en el plano.

## 🎯 Características

- **Catálogo de configuraciones predefinidas de 3 a 6 discos**
- **Análisis variacional**: Matriz de contacto, rolling space, gradiente proyectado, Hessiano intrínseco
- **Visualización**: Gráficos de discos, grafos de contacto, espectros
- **API**: Funciones de alto nivel para análisis rápido

## 📦 Instalación

### Desde PyPI (AUN NO FUNCIONA)
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

## 📖 Documentación

- **[Docs](docs/index.md)** - Guía exhaustiva, incluyendo fundamentos matemáticos
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
disk-packing-analysis/
├── extremal_packings/          # Código fuente principal
│   ├── __init__.py            # API pública
│   ├── analysis.py            # Pipeline de análisis
│   ├── catalog.py             # Catálogo de configuraciones
│   ├── configurations.py      # Clase Configuration
│   ├── constraints.py         # Matriz A y rolling space
│   ├── contact_graphs.py      # Validación de grafos
│   ├── hessian.py            # Hessiano K y H
│   ├── interface.py          # Funciones de alto nivel
│   ├── json_loader.py        # Carga desde JSON
│   ├── perimeter.py          # Perímetros y convex hull
│   └── plotting.py           # Visualización
├── data/                      # Configuraciones JSON
│   ├── 3disks.json
│   ├── 4disks.json
│   ├── 5disks.json
│   └── 6disks.json
├── tests/                     # Suite de tests
│   ├── test_analysis.py
│   ├── test_catalog.py
│   ├── test_configurations.py
│   ├── test_constraints.py
│   ├── test_contact_graphs.py
│   ├── test_hessian.py
│   └── test_perimeter.py
├── examples/                  # Ejemplos de uso
│   ├── basic_usage.py
│   └── advanced_usage.py
├── docs/                      # Documentación
│   ├── index.md
│   ├── api.md
│   └── DETAILED_DOCUMENTATION.md
├── pyproject.toml            # Configuración del proyecto
├── setup.py                  # Setup alternativo
├── README.md                 # Readme principal
└── LICENSE                   # Licencia MIT
```

## 🧪 Testing

```bash
pytest tests/
```

## 📝 Licencia

MIT License - Ver [LICENSE](LICENSE) para detalles.

## 👤 Autores

- **Fabián Andrés Henry Vilaxa**
- **Jose Ayala Hoffman**