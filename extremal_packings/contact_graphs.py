from __future__ import annotations
from typing import List, Tuple
from collections import deque


Edge = Tuple[int, int]


def is_connected(n: int, edges: List[Edge]) -> bool:
    if n == 0:
        return True
    if not edges:
        return n == 1  # si hay más de 1 vértice sin aristas, no es conexo

    adj = [[] for _ in range(n)]
    for (i, j) in edges:
        adj[i].append(j)
        adj[j].append(i)

    visited = [False] * n
    q = deque([0])
    visited[0] = True

    while q:
        u = q.popleft()
        for v in adj[u]:
            if not visited[v]:
                visited[v] = True
                q.append(v)

    return all(visited[i] or not adj[i] for i in range(n))


def max_degree(n: int, edges: List[Edge]) -> int:
    deg = [0] * n
    for (i, j) in edges:
        deg[i] += 1
        deg[j] += 1
    return max(deg) if n > 0 else 0


def check_graph_validity(n: int, edges: List[Edge], max_deg: int = 6) -> None:
    """
    Verifica condiciones básicas del grafo de contacto:
    - índices válidos
    - conectividad (salvo el caso trivial n=1)
    - grado máximo ≤ max_deg
    """
    for (i, j) in edges:
        if i < 0 or j < 0 or i >= n or j >= n:
            raise ValueError(f"Arista {(i, j)} tiene índice fuera de rango")

    if n > 1 and not is_connected(n, edges):
        raise ValueError("El grafo de contacto no es conexo")

    if max_degree(n, edges) > max_deg:
        raise ValueError(
            f"El grafo tiene vértices con grado > {max_deg}. "
            "No es compatible con contactos de discos unitarios."
        )