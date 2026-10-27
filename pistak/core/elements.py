"""Dati per elemento: numeri atomici, masse, raggi covalenti, raggi di van der Waals e colori di visualizzazione."""
from __future__ import annotations

# Simboli in ordine di numero atomico (fino al curio, come la tabella dei raggi covalenti).
SYMBOLS: tuple[str, ...] = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
)
ATOMIC_NUMBERS: dict[str, int] = {e: i + 1 for i, e in enumerate(SYMBOLS)}

# Masse atomiche medie in u (pesi atomici standard IUPAC, valori convenzionali; per gli elementi senza isotopi
# stabili il numero di massa dell'isotopo più longevo).
ATOMIC_MASSES: dict[str, float] = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Ne": 20.180, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974, "S": 32.06, "Cl": 35.45,
    "Ar": 39.95, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.922,
    "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224, "Nb": 92.906,
    "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82,
    "Sn": 118.71, "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91, "Ba": 137.33, "La": 138.91,
    "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97, "Hf": 178.49, "Ta": 180.95,
    "W": 183.84, "Re": 186.21, "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59, "Tl": 204.38,
    "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0, "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0,
    "Th": 232.04, "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0, "Cm": 247.0,
}

# Raggi covalenti in Å (B. Cordero et al., Dalton Trans. 2008, 2832). Per Mn, Fe e Co il valore a basso
# spin, per il carbonio quello sp³: con la tolleranza usata per i legami la differenza non conta.
COVALENT_RADII: dict[str, float] = {
    "H": 0.31, "He": 0.28, "Li": 1.28, "Be": 0.96, "B": 0.84, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "Ne": 0.58,
    "Na": 1.66, "Mg": 1.41, "Al": 1.21, "Si": 1.11, "P": 1.07, "S": 1.05, "Cl": 1.02, "Ar": 1.06,
    "K": 2.03, "Ca": 1.76, "Sc": 1.7, "Ti": 1.6, "V": 1.53, "Cr": 1.39, "Mn": 1.39, "Fe": 1.32, "Co": 1.26,
    "Ni": 1.24, "Cu": 1.32, "Zn": 1.22, "Ga": 1.22, "Ge": 1.2, "As": 1.19, "Se": 1.2, "Br": 1.2, "Kr": 1.16,
    "Rb": 2.2, "Sr": 1.95, "Y": 1.9, "Zr": 1.75, "Nb": 1.64, "Mo": 1.54, "Tc": 1.47, "Ru": 1.46, "Rh": 1.42,
    "Pd": 1.39, "Ag": 1.45, "Cd": 1.44, "In": 1.42, "Sn": 1.39, "Sb": 1.39, "Te": 1.38, "I": 1.39, "Xe": 1.4,
    "Cs": 2.44, "Ba": 2.15, "La": 2.07, "Ce": 2.04, "Pr": 2.03, "Nd": 2.01, "Pm": 1.99, "Sm": 1.98, "Eu": 1.98,
    "Gd": 1.96, "Tb": 1.94, "Dy": 1.92, "Ho": 1.92, "Er": 1.89, "Tm": 1.9, "Yb": 1.87, "Lu": 1.87, "Hf": 1.75,
    "Ta": 1.7, "W": 1.62, "Re": 1.51, "Os": 1.44, "Ir": 1.41, "Pt": 1.36, "Au": 1.36, "Hg": 1.32, "Tl": 1.45,
    "Pb": 1.46, "Bi": 1.48, "Po": 1.4, "At": 1.5, "Rn": 1.5, "Fr": 2.6, "Ra": 2.21, "Ac": 2.15, "Th": 2.06,
    "Pa": 2.0, "U": 1.96, "Np": 1.9, "Pu": 1.87, "Am": 1.8, "Cm": 1.69,
}

# Raggi di van der Waals in Å (A. Bondi, J. Phys. Chem. 1964, 68, 441; Mantina et al. 2009 per gli altri).
VDW_RADII: dict[str, float] = {
    "H": 1.10, "He": 1.40, "Li": 1.81, "Be": 1.53, "B": 1.92, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "Ne": 1.54,
    "Na": 2.27, "Mg": 1.73, "Al": 1.84, "Si": 2.10, "P": 1.80, "S": 1.80, "Cl": 1.75, "Ar": 1.88,
    "K": 2.75, "Ca": 2.31, "Ga": 1.87, "Ge": 2.11, "As": 1.85, "Se": 1.90, "Br": 1.85, "Kr": 2.02,
    "Rb": 3.03, "Sr": 2.49, "In": 1.93, "Sn": 2.17, "Sb": 2.06, "Te": 2.06, "I": 1.98, "Xe": 2.16,
    "Cs": 3.43, "Ba": 2.68, "Tl": 1.96, "Pb": 2.02, "Bi": 2.07,
    "Ni": 1.63, "Cu": 1.40, "Zn": 1.39, "Pd": 1.63, "Ag": 1.72, "Cd": 1.58, "Pt": 1.75, "Au": 1.66, "Hg": 1.55,
}

# Colori Jmol (RGB 0–255), per gli elementi più comuni; gli altri sono grigio-rosa.
COLORS: dict[str, tuple[int, int, int]] = {
    "H": (255, 255, 255), "B": (255, 181, 181), "C": (144, 144, 144), "N": (48, 80, 248), "O": (255, 13, 13),
    "F": (144, 224, 80), "Na": (171, 92, 242), "Mg": (138, 255, 0), "Al": (191, 166, 166), "Si": (240, 200, 160),
    "P": (255, 128, 0), "S": (255, 255, 48), "Cl": (31, 240, 31), "K": (143, 64, 212), "Ca": (61, 255, 0),
    "Mn": (156, 122, 199), "Fe": (224, 102, 51), "Co": (240, 144, 160), "Ni": (80, 208, 80), "Cu": (200, 128, 51),
    "Zn": (125, 128, 176), "Br": (166, 41, 41), "Ru": (36, 143, 143), "Rh": (10, 125, 140), "Pd": (0, 105, 133),
    "Ag": (192, 192, 192), "I": (148, 0, 148), "Pt": (208, 208, 224), "Au": (255, 209, 35), "Hg": (184, 184, 208),
    "Se": (255, 161, 0), "Sn": (102, 128, 128), "Ir": (23, 84, 135), "Re": (38, 125, 171), "Os": (38, 102, 150),
}
DEFAULT_COLOR = (255, 192, 203)
# Il carbonio della seconda molecola: ciano, come nella vista 3D del sito ("cyanCarbon").
CYAN_CARBON = (64, 208, 224)

DEFAULT_COVALENT = 1.5
DEFAULT_VDW = 1.9


def is_known(el: str) -> bool:
    return el in COVALENT_RADII


def atomic_number(el: str) -> int:
    """Numero atomico (0 se l'elemento non è nella tabella)."""
    return ATOMIC_NUMBERS.get(el, 0)


def atomic_mass(el: str) -> float:
    """Massa atomica in u (NaN se l'elemento non è nella tabella)."""
    return ATOMIC_MASSES.get(el, float("nan"))


def covalent_radius(el: str) -> float:
    return COVALENT_RADII.get(el, DEFAULT_COVALENT)


def vdw_radius(el: str) -> float:
    return VDW_RADII.get(el, DEFAULT_VDW)


def color(el: str) -> tuple[int, int, int]:
    return COLORS.get(el, DEFAULT_COLOR)
