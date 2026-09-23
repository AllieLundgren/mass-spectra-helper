"""
m/z reverse lookup over all_species.py — plain Python, no Streamlit.

Kept out of app.py on purpose: Python imports a module once per server process,
so the index and isotope patterns built here survive Streamlit reruns without
st.cache_* (which hung on Streamlit Cloud).

Only compositions and mass ranges are indexed up front (well under a second).
A search computes isotope patterns just for the species whose mass range can
reach the target m/z and remembers them for later searches.  Patterns are
convolutions of per-element patterns, grouped by nominal mass the same way
molmass' spectrum() does it.
"""

from __future__ import annotations

import re
import threading

from molmass import Formula
from molmass.elements import ELEMENTS

from all_species import all_species

PEAK_MIN_FRACTION   = 1e-5                            # peaks below 0.001 % are ignored
BACKGROUND_ELEMENTS = frozenset({"O", "H", "C", "N"})  # allowed with "Only these" when ticked

_SIMPLE_FORMULA = re.compile(r'([A-Z][a-z]?)(\d*)')
_index: dict | None = None
_index_lock = threading.Lock()
_element_cache: dict = {}  # (symbol, count) → (first_mass_number, abundances, abundance × mass)
_peak_cache: dict = {}     # composition key → [(mass_number, mz, abundance), ...]


def composition(formula_str: str) -> dict[str, int]:
    """Element counts for a formula, e.g. 'Na2Cl' → {'Na': 2, 'Cl': 1}. Raises on invalid input."""
    s = formula_str.strip()
    parts = _SIMPLE_FORMULA.findall(s)
    if parts and "".join(sym + n for sym, n in parts) == s:
        comp: dict[str, int] = {}
        for sym, n in parts:
            comp[sym] = comp.get(sym, 0) + (int(n) if n else 1)
    else:
        # Parentheses etc. (e.g. "Al(OH)3") — let molmass parse it
        comp = {sym: int(item.count) for sym, item in Formula(s).composition().items()}
    for sym in comp:
        ELEMENTS[sym]  # raises for non-elements such as the old "RE2O3"
    return comp


def _natural_isotopes(symbol: str) -> dict:
    isos = {mn: iso for mn, iso in ELEMENTS[symbol].isotopes.items() if iso.abundance > 0}
    if not isos:
        raise ValueError(f"{symbol} has no natural isotopes")
    return isos


def _combine(a, b):
    """Convolve two (first_mass_number, abundances, abundance × exact mass) patterns."""
    m0a, pa, wa = a
    m0b, pb, wb = b
    p = [0.0] * (len(pa) + len(pb) - 1)
    w = [0.0] * len(p)
    for i, (pi, wi) in enumerate(zip(pa, wa)):
        if pi == 0.0:
            continue
        for j, (pj, wj) in enumerate(zip(pb, wb)):
            p[i + j] += pi * pj
            w[i + j] += wi * pj + pi * wj
    return m0a + m0b, p, w


def _element_pattern(symbol: str, count: int):
    key = (symbol, count)
    if key in _element_cache:
        return _element_cache[key]
    if count == 1:
        isos = _natural_isotopes(symbol)
        m0 = min(isos)
        p = [0.0] * (max(isos) - m0 + 1)
        w = [0.0] * len(p)
        for mn, iso in isos.items():
            p[mn - m0] += iso.abundance
            w[mn - m0] += iso.abundance * iso.mass
        total = sum(p)
        out = (m0, [x / total for x in p], [x / total for x in w])
    else:
        half = _element_pattern(symbol, count // 2)
        out = _combine(half, half)
        if count % 2:
            out = _combine(out, _element_pattern(symbol, 1))
    _element_cache[key] = out
    return out


def species_peaks(comp_key: tuple) -> list[tuple[int, float, float]]:
    """[(mass_number, mz, abundance), ...] for a sorted composition tuple, neutral masses."""
    peaks = _peak_cache.get(comp_key)
    if peaks is None:
        out = None
        for sym, n in comp_key:
            pat = _element_pattern(sym, n)
            out = pat if out is None else _combine(out, pat)
        m0, p, w = out
        peaks = [(m0 + k, w[k] / pk, pk) for k, pk in enumerate(p) if pk >= PEAK_MIN_FRACTION]
        _peak_cache[comp_key] = peaks
    return peaks


def _build_index() -> dict:
    seen: set = set()
    invalid: set = set()
    mass_range: dict[str, tuple[float, float]] = {}
    index = {
        "names": [], "group": [], "keys": [], "elements": [],
        "lo": [], "hi": [], "atoms": [],
    }
    for group, formulas in all_species.items():
        for formula in formulas:
            try:
                comp = composition(formula)
                key = tuple(sorted(comp.items()))
                if key in seen:
                    continue
                for sym in comp:
                    if sym not in mass_range:
                        masses = [iso.mass for iso in _natural_isotopes(sym).values()]
                        mass_range[sym] = (min(masses), max(masses))
            except Exception:
                invalid.add(formula)
                continue
            seen.add(key)
            index["names"].append(formula)
            index["group"].append(group)
            index["keys"].append(key)
            index["elements"].append(frozenset(comp))
            index["lo"].append(sum(mass_range[s][0] * n for s, n in key))
            index["hi"].append(sum(mass_range[s][1] * n for s, n in key))
            index["atoms"].append(sum(comp.values()))
    index["invalid"] = sorted(invalid)
    return index


def get_species_index() -> dict:
    """The species index, built on first use and shared by every session."""
    global _index
    if _index is None:
        with _index_lock:
            if _index is None:
                _index = _build_index()
    return _index


def _elements_match(els: frozenset, requested: set[str], mode: str, allow_background: bool) -> bool:
    if mode == "All of these":
        return requested <= els
    if mode == "Only these":
        allowed = (requested | BACKGROUND_ELEMENTS) if allow_background else requested
        return els <= allowed
    return bool(requested & els)  # "Any of these"


def lookup_peaks(
    target: float,
    tol: float,
    groups: set[str],
    charges: tuple[int, ...] = (1, 2),
    min_abundance: float = 0.0,
    elements: set[str] | None = None,
    mode: str = "Any of these",
    allow_background: bool = True,
) -> tuple[list[dict], int]:
    """
    Isotope peaks within ±tol of target that pass the group/charge/abundance/element filters.
    2+ is only considered for atoms and diatomics.
    Returns (match rows sorted by error, number of candidate species checked).
    """
    index = get_species_index()
    names, grp, keys, els = index["names"], index["group"], index["keys"], index["elements"]
    lo, hi, atoms = index["lo"], index["hi"], index["atoms"]
    rows: list[dict] = []
    n_checked = 0

    for z in charges:
        low, high = (target - tol) * z, (target + tol) * z  # compare in neutral-mass units
        for i in range(len(names)):
            if lo[i] > high or hi[i] < low or grp[i] not in groups:
                continue
            if z > 1 and atoms[i] > 2:
                continue
            if elements and not _elements_match(els[i], elements, mode, allow_background):
                continue
            n_checked += 1
            for mass_number, mz, abund in species_peaks(keys[i]):
                mz_z = mz / z
                err = abs(mz_z - target)
                if err <= tol and abund >= min_abundance:
                    rows.append({
                        "Species":       names[i],
                        "Isotope":       f"{names[i]}-{mass_number}" + (" (2+)" if z == 2 else ""),
                        "Charge":        z,
                        "m/z":           mz_z,
                        "Abundance (%)": abund * 100,
                        "Error":         err,
                        "Group":         grp[i].replace("_", " "),
                    })

    rows.sort(key=lambda r: r["Error"])
    return rows, n_checked
