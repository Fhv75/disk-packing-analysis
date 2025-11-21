from __future__ import annotations
from typing import List, Tuple
import json
import numpy as np
from .configurations import Configuration


def load_from_json(path: str) -> Configuration:
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    coords = np.array(data["coords"], dtype=float)
    edges: List[Tuple[int, int]] = [tuple(e) for e in data["edges"]]
    name = data.get("name")
    hull = data.get("hull_vertices")
    return Configuration(coords=coords, edges=edges, name=name, hull_vertices=hull)


def save_to_json(config: Configuration, path: str) -> None:
    data = {
        "name": config.name,
        "coords": config.coords.tolist(),
        "edges": config.edges,
        "hull_vertices": config.hull_vertices,
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)