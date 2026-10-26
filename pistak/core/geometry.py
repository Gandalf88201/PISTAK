"""Geometria elementare: piano medio, involucro convesso 2D, ritaglio di poligoni, angoli."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / (n or 1.0)


def angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    c = np.dot(u, v) / ((np.linalg.norm(u) * np.linalg.norm(v)) or 1.0)
    return float(np.degrees(np.arccos(np.clip(c, -1.0, 1.0))))


def acute_angle_deg(u: np.ndarray, v: np.ndarray) -> float:
    """Angolo tra due rette (0–90°): il verso dei vettori non conta."""
    a = angle_deg(u, v)
    return a if a <= 90 else 180 - a


@dataclass(frozen=True)
class Plane:
    centroid: np.ndarray
    normal: np.ndarray
    long: np.ndarray  # direzione di massima estensione nel piano
    rms: float  # scarto quadratico medio dei punti dal piano (Å)


def _canonical(v: np.ndarray) -> np.ndarray:
    """Segno fissato (componente maggiore positiva): risultati riproducibili tra un'esecuzione e l'altra."""
    return v if v[np.argmax(np.abs(v))] >= 0 else -v


def fit_plane(ps: np.ndarray) -> Plane:
    """Piano medio di un insieme di punti: normale = direzione di minima dispersione, asse lungo = di massima."""
    ps = np.asarray(ps, dtype=float)
    c = ps.mean(axis=0)
    d = ps - c
    values, vectors = np.linalg.eigh(d.T @ d)  # autovalori in ordine crescente
    return Plane(c, _canonical(vectors[:, 0]), _canonical(vectors[:, 2]), float(np.sqrt(max(0.0, values[0]) / len(ps))))


def convex_hull(points: np.ndarray) -> np.ndarray:
    """Involucro convesso di punti 2D (catena monotona di Andrew), in senso antiorario."""
    ps = sorted(map(tuple, np.asarray(points, dtype=float)))
    if len(ps) < 3:
        return np.array(ps)

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: list[tuple] = []
    for p in ps:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper: list[tuple] = []
    for p in reversed(ps):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return np.array(lower[:-1] + upper[:-1])


def points_inside(poly: np.ndarray, pts: np.ndarray) -> np.ndarray:
    """Quali punti 2D cadono dentro (o sul bordo, con tolleranza) di un poligono convesso antiorario."""
    pts = np.asarray(pts, dtype=float)
    if len(poly) < 3:
        return np.zeros(len(pts), dtype=bool)
    a = poly
    b = np.roll(poly, -1, axis=0)
    ab = b - a  # E×2
    ap = pts[:, None, :] - a[None, :, :]  # P×E×2
    cr = ab[None, :, 0] * ap[:, :, 1] - ab[None, :, 1] * ap[:, :, 0]
    return (cr >= -1e-9).all(axis=1)


def polygon_area(poly: np.ndarray) -> float:
    if len(poly) < 3:
        return 0.0
    x, y = poly[:, 0], poly[:, 1]
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2)


def clip_convex(subject: np.ndarray, clip: np.ndarray) -> np.ndarray:
    """Intersezione di due poligoni convessi antiorari (Sutherland–Hodgman). Vuota se non si toccano."""
    out = [tuple(p) for p in subject]
    for i in range(len(clip)):
        a, b = clip[i], clip[(i + 1) % len(clip)]
        if not out:
            break
        inp, out = out, []
        side = lambda p: (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0])  # noqa: E731, B023 (usata subito, nello stesso giro)
        for j in range(len(inp)):
            p, q = inp[j], inp[(j + 1) % len(inp)]
            sp, sq = side(p), side(q)
            if sp >= 0:
                out.append(p)
            if (sp >= 0) != (sq >= 0):
                t = sp / (sp - sq)
                out.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    return np.array(out)
