# Extremal Packings - Análisis de Empaquetamientos de Discos Congruentes

**Versión:** 1.0.0  
**Autor:** Fabián Andrés Henry Vilaxa, Jose Ayala Hoffman  
**Licencia:** MIT

---

## Tabla de Contenidos

1. [Introducción](#introducción)
2. [Fundamentos Matemáticos](#fundamentos-matemáticos)
3. [Arquitectura del Paquete](#arquitectura-del-paquete)
4. [Módulos Detallados](#módulos-detallados)
5. [Catálogo de Configuraciones](#catálogo-de-configuraciones)
6. [Ejemplos Avanzados](#ejemplos-avanzados)
7. [Guía de Desarrollo](#guía-de-desarrollo)
8. [Referencias](#referencias)

---

## Introducción

### Motivación

El paquete `extremal_packings` surge de la investigación sobre **perímetros extremos en empaquetamientos de discos congruentes**. El objetivo es estudiar configuraciones de discos unitarios tangentes que minimizan o maximizan el perímetro de su envolvente convexa.

### Alcance

El paquete proporciona:

- **Análisis geométrico**: Envolventes convexas, perímetros, grafos de contacto
- **Análisis variacional**: Rolling space, Análisis de criticidad, Análisis de estabilidad.
- **Catálogo**: Configuraciones predefinidas de 3 a 6 discos
- **Visualización**: Gráficos estáticos con Matplotlib

### Instalación

```bash
# Desde PyPI
pip install extremal-packings

# Desde el repositorio
git clone https://github.com/fhenr/disk-packing-analysis.git
cd disk-packing-analysis
pip install -e .
```

### Dependencias

- **Python**: ≥3.8
- **NumPy**: ≥1.20.0 (álgebra lineal)
- **SciPy**: ≥1.7.0 (convex hull, SVD)
- **Matplotlib**: ≥3.3.0 (visualización)
- **NetworkX**: ≥2.5 (grafos de contacto)

---

## Fundamentos Matemáticos

### 1. Configuración de Discos

**Definición**: Una configuración $\mathcal{C}$ consiste en:

- $n$ discos unitarios (radio $r = 1$)
- Centros $c_1, \ldots, c_n \in \mathbb{R}^2$
- Grafo de contacto $G = (V, E)$ donde $(i,j) \in E \iff \|c_j - c_i\| = 2$

**Espacio de configuraciones**: $\mathcal{D} = \mathbb{R}^{2n}$ (coordenadas apiladas: $[c_1^x, c_1^y, \ldots, c_n^x, c_n^y]$)

### 2. Matriz de Contacto $A(c)$

Para cada contacto $(i,j) \in E$, definimos el vector unitario:

$$u_{ij} = \frac{c_j - c_i}{\|c_j - c_i\|} = \frac{c_j - c_i}{2}$$

La matriz de contacto $A(c) \in \mathbb{R}^{m \times 2n}$ tiene una fila por contacto:

$$\text{fila}_k = [\ldots, 0, -u_{ij}, 0, \ldots, 0, u_{ij}, 0, \ldots]$$

con $-u_{ij}$ en las posiciones $(2i, 2i+1)$ y $u_{ij}$ en $(2j, 2j+1)$.

**Propiedad**: $A(c) \cdot \delta{c} = 0 \iff$ la deformación $\delta{c}$ preserva todos los contactos infinitesimalmente.

### 3. Rolling Space

El **rolling space** es el kernel de $A(c)$:

$$\text{Roll}(c) = \ker(A(c)) = \{\delta{c} \in \mathbb{R}^{2n} : A(c) \cdot \delta{c} = 0\}$$

**Dimensión**:

$$\dim(\text{Roll}(c)) = 2n - \text{rank}(A(c))$$

**Interpretación física**:
- Deformaciones infinitesimales que mantienen las distancias de contacto constantes.
- Incluye movimientos rígidos: 2 traslaciones + 1 rotación = 3 grados.
- Si $\dim(\text{Roll}) = 3$, la configuración es **infinitesimalmente rígida**.
- Si $\dim(\text{Roll}) > 3$, hay **flexibilidad** (deformaciones no triviales).

### 4. Funcional de Perímetro

El perímetro del cluster de discos es:

$$P(c) = \text{perímetro del casco convexo de } \bigcup_{i=1}^n \{x : \|x - c_i\| \leq 1\}$$

**Simplificación**: Para configuraciones con discos en la frontera bien definida, se puede calcular mediante:

1. Casco convexo de centros: $\text{CH}(c_1, \ldots, c_n)$
2. Agregar arcos circulares en los vértices del casco

### 5. Hessiano del Perímetro

#### 5.1 Hessiano Global $K(c)$

Matriz simétrica $K(c) \in \mathbb{R}^{2n \times 2n}$ que representa la segunda derivada de $P(c)$:

$$K(c) = \nabla^2 P(c)$$

**Estructura por bloques**:

$$K(c) = \begin{pmatrix}
K_{11} & \cdots & K_{1n} \\
\vdots & \ddots & \vdots \\
K_{n1} & \cdots & K_{nn}
\end{pmatrix}$$

donde cada $K_{ij} \in \mathbb{R}^{2 \times 2}$.

**Contribuciones**:
- **Diagonales** $K_{ii}$: Interacción del disco $i$ con todos sus vecinos
- **Fuera de diagonal** $K_{ij}$ ($i \neq j$): Si $(i,j) \in E$, hay contribución del contacto

#### 5.2 Hessiano Intrínseco $H$

Proyección del Hessiano global al rolling space:

$$H = R^T K(c) R$$

donde $R \in \mathbb{R}^{2n \times d}$ es una base ortonormal de $\text{Roll}(c)$ con $d = \dim(\text{Roll})$.

**Propiedades**:
- $H \in \mathbb{R}^{d \times d}$ es simétrica
- Autovalores $\lambda_1 \leq \lambda_2 \leq \cdots \leq \lambda_d$

**Interpretación**:
- $\lambda_i > 0$: Curvatura positiva en dirección $i$ (mínimo local estable)
- $\lambda_i = 0$: Dirección neutra (degeneración)
- $\lambda_i < 0$: Curvatura negativa (punto de silla, no es mínimo)

### 6. Condiciones de Optimalidad

Para que $c$ sea un punto crítico del perímetro en el rolling space:

1. **Primera variación**: $\nabla P(c) \perp \text{Roll}(c)$ (gradiente ortogonal al rolling space). Equivalente a $\langle \nabla \text{Per}(c), \; \delta{c} \rangle = 0 \quad \forall \; \delta{c} \in \text{Roll(c)}$.
2. **Segunda variación**: $H = R^T K(c) R$ debe ser semidefinida positiva ($\lambda_i \geq 0 \, \forall i$)

---

## Arquitectura del Paquete

### Estructura de Directorios

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

### Flujo de Datos

```
Configuration (coords, edges)
    ↓
[Validación]
    ↓
A(c) ← build_contact_matrix
    ↓
R ← rolling_space_basis(A)
    ↓
K(c) ← build_unconstrained_hessian
    ↓
H = R^T K R ← project_to_roll
    ↓
λ₁, ..., λ_d ← intrinsic_spectrum(H)
    ↓
AnalysisResult
```

---

## Módulos Detallados

### 1. `configurations.py`

Define la clase base `Configuration`.

#### Clase `Configuration`

```python
@dataclass
class Configuration:
    coords: np.ndarray          # (n, 2) coordenadas de centros
    edges: list[tuple[int, int]] # Lista de contactos
    name: str = ""              # Nombre identificador
    radius: float = 1.0         # Radio (siempre 1.0)
```

**Métodos**:

- `n` (property): Número de discos
- `__post_init__()`: Validación de datos de entrada

**Ejemplo**:

```python
coords = np.array([[0, 0], [2, 0], [1, np.sqrt(3)]])
edges = [(0, 1), (1, 2), (2, 0)]
config = Configuration(coords=coords, edges=edges, name="Triangle")
```

---

### 2. `contact_graphs.py`

Validación de grafos de contacto.

#### Función `check_graph_validity`

```python
def check_graph_validity(n: int, edges: list[tuple[int, int]]) -> None:
    """Valida que el grafo de contacto sea físicamente realizable."""
```

**Validaciones**:

1. **Conectividad**: El grafo debe ser conexo (componente única)
2. **Grados**: Ningún vértice puede tener grado > 6 (*kissing number* en $\mathbb{R}^2$)
3. **Índices**: Todos los índices deben estar en [0, n-1]

**Excepciones**:
- `ValueError` si el grafo no cumple alguna condición

**Ejemplo**:

```python
check_graph_validity(4, [(0,1), (1,2), (2,3), (3,0)])  # OK
check_graph_validity(4, [(0,1), (2,3)])  # ValueError: desconectado
```

---

### 3. `constraints.py`

Construcción de la matriz de contacto y rolling space.

#### Función `build_contact_matrix`

```python
def build_contact_matrix(config: Configuration) -> np.ndarray:
    """Construye matriz A(c) de dimensión (m, 2n)."""
```

**Algoritmo**:

```python
Para cada contacto (i, j):
    1. u_ij = (c_j - c_i) / ||c_j - c_i||
    2. fila_k = [0, ..., -u_ij[0], -u_ij[1], ..., u_ij[0], u_ij[1], ..., 0]
```

**Salida**: Matriz $A \in \mathbb{R}^{m \times 2n}$ con filas ortonormales (si contactos son ortogonales).

#### Función `rolling_space_basis`

```python
def rolling_space_basis(A: np.ndarray) -> np.ndarray:
    """Calcula base ortonormal del ker(A) mediante SVD."""
```

**Algoritmo (SVD)**:

1. Descomposición: $A = U \Sigma V^T$
2. Valores singulares: $\sigma_1 \geq \cdots \geq \sigma_r > 0$ (rank = r)
3. Rolling space: Columnas de $V$ correspondientes a $\sigma_i < \text{tol}$
4. Dimensión: $d = 2n - r$

**Tolerancia**: $\text{tol} = 10^{-10}$ (valores singulares menores se consideran 0)

**Ejemplo**:

```python
A = build_contact_matrix(config)
R = rolling_space_basis(A)  # R: (2n, d)
print(f"Dimension rolling space: {R.shape[1]}")
```

---

### 4. `hessian.py`

Cálculo del Hessiano global y proyección al rolling space.

#### Función `build_unconstrained_hessian`

```python
def build_unconstrained_hessian(config: Configuration) -> np.ndarray:
    """Construye Hessiano global K(c) en R^{2n × 2n}."""
```

**Estructura**:

$$K = \begin{pmatrix}
K_{11} & \cdots & K_{1n} \\
\vdots & \ddots & \vdots \\
K_{n1} & \cdots & K_{nn}
\end{pmatrix}, \quad K_{ij} \in \mathbb{R}^{2 \times 2}$$

**Contribuciones por contacto** $(i, j)$:

$$K_{ij} = -\frac{1}{2} P_{ij}, \quad K_{ii} += \frac{1}{2} P_{ij}, \quad K_{jj} += \frac{1}{2} P_{ij}$$

donde $P_{ij} = I - u_{ij} u_{ij}^T$ es el proyector perpendicular a $u_{ij}$.

#### Función `project_to_roll`

```python
def project_to_roll(K: np.ndarray, R: np.ndarray) -> np.ndarray:
    """Proyecta K al rolling space: H = R^T K R."""
```

**Operación**: Multiplicación matricial $H = R^T K R$

**Dimensiones**:
- Entrada: $K \in \mathbb{R}^{2n \times 2n}$, $R \in \mathbb{R}^{2n \times d}$
- Salida: $H \in \mathbb{R}^{d \times d}$

#### Función `intrinsic_spectrum`

```python
def intrinsic_spectrum(H: np.ndarray) -> np.ndarray:
    """Calcula autovalores de H ordenados de menor a mayor."""
```

**Método**: `numpy.linalg.eigvalsh` (para matrices simétricas)

**Salida**: Array 1D con autovalores $\lambda_1 \leq \cdots \leq \lambda_d$

---

### 5. `perimeter.py`

Cálculo de perímetros y envolventes convexas.

#### Función `compute_hull`

```python
def compute_hull(config: Configuration) -> np.ndarray:
    """Calcula casco convexo de centros usando scipy.spatial.ConvexHull."""
```

**Salida**: Índices de vértices en el casco, ordenados en sentido antihorario.

#### Función `perimeter_centers`

```python
def perimeter_centers(config: Configuration) -> float:
    """Perímetro del casco convexo de centros."""
```

**Algoritmo**:

1. Calcular casco: `hull = ConvexHull(coords)`
2. Sumar longitudes: $\sum_{k} \|v_{k+1} - v_k\|$

#### Función `perimeter_disks`

```python
def perimeter_disks(config: Configuration) -> float:
    """Perímetro del cluster de discos (centros + arcos circulares)."""
```

**Algoritmo**:

1. Casco de centros con índices $h_1, \ldots, h_k$
2. Para cada arista $(h_i, h_{i+1})$:
   - Segmento recto si $(h_i, h_{i+1}) \in E$
   - Arco circular de ángulo $\theta$ si no hay contacto
3. Perímetro total = suma de segmentos + arcos

**Complejidad**: $O(n \log n)$ (dominado por ConvexHull)

---

### 6. `analysis.py`

Pipeline completo de análisis.

#### Clase `AnalysisResult`

Contenedor de resultados con propiedades calculadas:

```python
@property
def rolling_dim(self) -> int:
    """Dimensión del rolling space."""
    return self.R.shape[1]

@property
def is_rigid(self) -> bool:
    """True si dim(Roll) ≤ 3 (solo movimientos rígidos)."""
    return self.rolling_dim <= 3

@property
def is_flexible(self) -> bool:
    """True si dim(Roll) > 3 (hay flexibilidad)."""
    return not self.is_rigid

@property
def has_negative_eigenvalue(self) -> bool:
    """True si existe λ < -10^{-10}."""
    if len(self.eigenvalues) == 0:
        return False
    return self.eigenvalues[0] < -1e-10
```

#### Función `analyze_configuration`

```python
def analyze_configuration(config: Configuration) -> AnalysisResult:
    """Pipeline completo: validación → A → R → K → H → λ."""
```

**Pasos**:

1. Validación del grafo
2. Matriz de contacto $A(c)$
3. Rolling space $R = \ker(A)$
4. Perímetros (centros y discos)
5. Gradiente $\nabla{\text{Per}}(c)$ y proyección sobre $\text{Roll}(c)$
6. Hessiano global $K(c)$
7. Hessiano intrínseco $H = R^T K R$
8. Espectro de $H$

**Ejemplo completo**:

```python
from extremal_packings import load_configuration, analyze_configuration

config = load_configuration("D5-7")
result = analyze_configuration(config)

print(f"n = {config.n}, m = {len(config.edges)}")
print(f"dim(Roll) = {result.rolling_dim}")
print(f"Rígida: {result.is_rigid}")
print(f"Crítica: {result.is_critical}")
print(f"Autovalores: {result.eigenvalues}")
print(f"Perímetro: {result.perimeter_disks:.4f}")
```

---

### 7. `catalog.py`

Gestión del catálogo de configuraciones predefinidas.

#### Sistema de Nomenclatura

Formato: **`D{n}-{idx}`**

- `n`: Número de discos (3-6)
- `idx`: Índice secuencial (1, 2, 3, ...)

Ejemplos:
- `D3-1`: Triángulo equilátero
- `D3-2`: Cadena de 3 discos
- `D5-7`: Pentágono regular
- `D6-46`: Configuración colineal de 6 discos

#### Funciones Principales

```python
def list_configurations() -> list[str]:
    """Lista todos los nombres disponibles, ordenados naturalmente."""

def load_configuration(name: str) -> Configuration:
    """Carga configuración por nombre."""

def get_configurations_by_size(n: int) -> list[str]:
    """Filtra configuraciones con n discos."""

def get_catalog_stats() -> dict[str, int]:
    """Estadísticas: total, por tamaño, rango."""
```

---

### 8. `json_loader.py`

Carga de configuraciones desde archivos JSON.

#### Formato JSON

```json
{
    "version": "1.1",
    "indexing": "0-based",
    "angles": "degrees",
    "radius": "1",
    "graphs": [
        {
            "discos": 3,
            "centros": [[0,0], [2,0], ["1","sqrt(3)"]],
            "contactos": [[0,1], [1,2], [2,0]]
        }
    ]
}
```

**Características**:

- Soporte para expresiones simbólicas: `"sqrt(3)"`, `"2*cosd(60)"`
- Evaluación con `sympy`
- Validación de índices de contactos

#### Función `load_all_configurations`

```python
def load_all_configurations(data_dir: str) -> dict[str, Configuration]:
    """Carga todas las configs de data/*.json."""
```

**Proceso**:

1. Escanear `*.json` en `data_dir`
2. Para cada archivo, parsear JSON
3. Evaluar expresiones simbólicas
4. Crear objetos `Configuration`
5. Asignar nombres: `D{n}-{idx}`

---

### 9. `plotting.py`

Visualización con Matplotlib.

#### Función `plot_disks`

```python
def plot_disks(
    config: Configuration,
    show_hull: bool = True,
    ax: Optional[plt.Axes] = None
) -> plt.Figure:
    """Grafica discos con envolvente convexa opcional."""
```

**Elementos visualizados**:

- Círculos: Discos con borde negro
- Centros: Puntos rojos
- Contactos: Líneas punteadas grises
- Casco convexo: Polígono azul semitransparente (si `show_hull=True`)

#### Función `plot_contact_graph`

```python
def plot_contact_graph(
    config: Configuration,
    show_normals: bool = False,
    ax: Optional[plt.Axes] = None
) -> plt.Figure:
    """Visualiza grafo de contacto con normales opcionales."""
```

**Elementos**:

- Nodos: Círculos en posiciones de centros
- Aristas: Líneas entre contactos
- Normales: Vectores unitarios $u_{ij}$ (si `show_normals=True`)

#### Función `plot_spectrum`

```python
def plot_spectrum(
    eigenvalues: np.ndarray,
    title: str = "",
    ax: Optional[plt.Axes] = None
) -> plt.Figure:
    """Gráfico de barras de autovalores."""
```

**Formato**:

- Eje X: Índice del autovalor
- Eje Y: Valor de $\lambda_i$
- Colores: Rojo (< 0), Amarillo (≈ 0), Verde (> 0)

---

### 10. `interface.py`

Funciones de alto nivel para usuarios finales.

#### Función `print_analysis_summary`

```python
def print_analysis_summary(result: AnalysisResult) -> None:
    """Imprime resumen formateado del análisis."""
```

**Salida ejemplo**:

```
========================================
Configuración: D5-7
========================================
Discos: 5
Contactos: 5
Dimensión rolling space: 5
Rígida: False

Autovalores del Hessiano intrínseco:
  λ[0] =  0.0000e+00
  λ[1] =  0.0000e+00
  λ[2] =  6.6910e-01
  λ[3] =  1.0000e+00
  λ[4] =  1.0000e+00

Perímetros:
  Centros: 10.0000
  Discos:  16.2832
========================================
```

#### Función `show_complete_analysis`

```python
def show_complete_analysis(config: Configuration) -> None:
    """Análisis completo con 3 gráficos en una figura."""
```

**Layout**: 1 fila × 3 columnas

1. **Izquierda**: Discos + casco (`plot_disks`)
2. **Centro**: Grafo de contacto (`plot_contact_graph`)
3. **Derecha**: Espectro (`plot_spectrum`)

#### Función `create_dashboard`

```python
def create_dashboard(configs: list[str]) -> None:
    """Dashboard comparativo para múltiples configuraciones."""
```

**Tabla comparativa**:

| Config | n | m | dim(Roll) | Perímetro | λ_min | λ_max |
|--------|---|---|-----------|-----------|-------|-------|
| D5-1   | 5 | 4 | 5         | 14.5678   | 0 | 0.83  |
| D5-7   | 5 | 5 | 3         | 16.2832   | 0.00  | 1  |

---

## Catálogo de Configuraciones

### Estadísticas del Catálogo

```python
from extremal_packings import get_catalog_stats

stats = get_catalog_stats()
# stats = {
#     'total': 65,
#     'by_size': {3: 2, 4: 5, 5: 13, 6: 45},
#     'min_disks': 3,
#     'max_disks': 6
# }
```

### Configuraciones Destacadas

#### D3-1: Triángulo Equilátero

```python
config = load_configuration("D3-1")
# 3 discos, 3 contactos
# Rígida (dim(Roll) = 3)
# Perímetro mínimo para 3 discos
```

#### D3-2: Cadena de 3 Discos

```python
config = load_configuration("D3-2")
# 3 discos, 2 contactos
# Flexible (dim(Roll) = 4)
```

#### D5-7: Pentágono Regular

```python
config = load_configuration("D5-7")
# 5 discos en vértices de pentágono regular
# 5 contactos formando ciclo
# Rígida pero con λ = 0 (frontera de estabilidad)
```

#### D6-1 a D6-45: Hexágono y Variantes

- 46 configuraciones distintas de 6 discos

---

## Ejemplos Avanzados

### Ejemplo 1: Análisis Comparativo

```python
from extremal_packings import *

configs = ["D5-1", "D5-7", "D5-11"]

print(f"{'Config':<10} {'Roll dim':<10} {'Perímetro':<12} {'Mín λ':<12}")
print("-" * 50)

for name in configs:
    config = load_configuration(name)
    result = analyze_configuration(config)
    
    print(f"{name:<10} {result.rolling_dim:<10} "
          f"{result.perimeter_disks:<12.4f} "
          f"{result.eigenvalues[0]:<12.4e}")
```



### Ejemplo 2: Pipeline Personalizado

```python
from extremal_packings.constraints import build_contact_matrix, rolling_space_basis
from extremal_packings.hessian import build_unconstrained_hessian, project_to_roll

# Crear configuración
config = Configuration(
    coords=np.array([[0,0], [2,0], [4,0]]),
    edges=[(0,1), (1,2)]
)

# Pipeline manual
A = build_contact_matrix(config)
R = rolling_space_basis(A)
K = build_unconstrained_hessian(config)
H = project_to_roll(K, R)

print(f"A: {A.shape}, R: {R.shape}, K: {K.shape}, H: {H.shape}")
```

---

## Guía de Desarrollo

### Estructura de Tests

Tests organizados por módulo:

```
tests/
├── test_configurations.py   # Configuration class
├── test_contact_graphs.py   # Validación de grafos
├── test_constraints.py      # A(c) y rolling space
├── test_hessian.py         # K, H, espectro
├── test_perimeter.py       # Perímetros y convex hull
├── test_catalog.py         # Catálogo y carga
└── test_analysis.py        # Pipeline completo
```

### Ejecutar Tests

```bash
# Todos los tests
pytest tests/

# Con cobertura
pytest --cov=extremal_packings tests/

# Test específico
pytest tests/test_analysis.py::TestAnalyzeConfiguration::test_triangle
```

### Agregar Nueva Configuración

1. **Crear entrada JSON en `data/{n}disks.json`**:

```json
{
    "discos": 5,
    "centros": [[0,0], [2,0], ...],
    "contactos": [[0,1], [1,2], ...]
}
```

2. **Validar con test**:

```python
def test_new_config():
    config = load_configuration("D5-14")
    assert config.n == 5
    assert len(config.edges) > 0
```

3. **Analizar**:

```python
result = analyze_configuration(config)
print_analysis_summary(result)
```

---

## Referencias

### Artículos Científicos

1. **Connelly, R.** (1980). *The rigidity of polyhedral surfaces*. Mathematics Magazine, 52(5), 275-283.

2. **Thurston, W.** (1998). *Shapes of polyhedra and triangulations of the sphere*. The Epstein Birthday Schrift, 511-549.

3. **Gromov, M.** (1983). *Filling Riemannian manifolds*. Journal of Differential Geometry, 18(1), 1-147.

### Libros de Referencia

4. **Grünbaum, B., & Shephard, G. C.** (1987). *Tilings and Patterns*. W. H. Freeman.

5. **Pach, J., & Agarwal, P. K.** (1995). *Combinatorial Geometry*. Wiley-Interscience.

### Recursos Online

- **OEIS A000055**: Número de grafos conexos con n vértices
- **ConvexHull Documentation**: [SciPy Spatial](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html)

---

## Apéndices

### A. Nomenclatura Matemática

| Símbolo | Significado |
|---------|-------------|
| $n$ | Número de discos |
| $m$ | Número de contactos |
| $c_i$ | Centro del disco $i$ |
| $r$ | Radio (siempre 1) |
| $G = (V, E)$ | Grafo de contacto |
| $A(c)$ | Matriz de contacto ($m \times 2n$) |
| $\text{Roll}(c)$ | Rolling space ($\dim = 2n - \text{rank}(A)$) |
| $R$ | Base ortonormal de $\text{Roll}(c)$ |
| $K(c)$ | Hessiano global ($2n \times 2n$) |
| $H$ | Hessiano intrínseco ($d \times d$) |
| $\lambda_i$ | Autovalor $i$ de $H$ |
| $P(c)$ | Perímetro del cluster |

### B. Glosario de Términos

**Punto Crítico del Perímetro**: $\langle \nabla \text{Per}(c), \; \delta{c} \rangle = 0$ para todo $\delta{c} \in \text{Roll(c)}$

**Configuración rígida**: $\dim(\text{Roll}) \leq 3$ (solo movimientos rígidos).

**Configuración flexible**: $\dim(\text{Roll}) > 3$ (deformaciones no triviales).

**Mínimo local estable**: $H$ semidefinida positiva ($\lambda_i \geq 0 \, \forall i$).

**Punto de silla**: Existe $\lambda_i < 0$.

**Grafo de contacto**: Grafo $G = (V, E)$ donde $V = \{1, \ldots, n\}$ y $(i,j) \in E \iff \|c_j - c_i\| = 2$.

**Convex hull**: Polígono convexo más pequeño que contiene todos los centros.

**Rolling space**: Espacio tangente de configuraciones con contactos preservados.

---

**Fin del documento** | Versión 1.0.1 | Última actualización: 2025
