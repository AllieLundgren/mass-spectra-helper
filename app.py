"""
Exum Instruments — Isotope Plots
===============================================================

Run:
    streamlit run isotope_app.py
"""

from __future__ import annotations

import io
import re
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from molmass import Formula
from all_species import all_species
from collections import defaultdict


def _is_single_isotope(formula_str: str) -> bool:
    """
    Returns True if the input is a specific single-isotope (e.g. 57Fe, Mo-92, Fe-56),
    which would only produce one peak. We want families only.
    """
    s = formula_str.strip()
    # Leading-mass single element: "57Fe", "92Mo"
    if re.match(r'^\d+[A-Z][a-z]?$', s):
        return True
    # Dash notation single element: "Fe-56", "Mo-92"
    if re.match(r'^[A-Z][a-z]?-\d+$', s):
        return True
    # Bracket notation: "[56Fe]"
    if re.match(r'^\[\d+[A-Z][a-z]?\]$', s):
        return True
    return False

# ── Optional H5 support ───────────────────────────────────────────────────────
try:
    import h5py
    _H5_AVAILABLE = True
except ImportError:
    _H5_AVAILABLE = False

# ── Optional baseline correction ──────────────────────────────────────────────
try:
    from pybaselines import Baseline
    _BASELINE_AVAILABLE = True
except ImportError:
    _BASELINE_AVAILABLE = False


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Isotope pattern generation  (inlined from isotopes.py)
# ═════════════════════════════════════════════════════════════════════════════

def _parse_formula(formula_str: str) -> Formula:
    """
    Parse a formula string into a molmass Formula object.
    Accepts standard formulas (Fe, Fe2O3) and isotope-specific
    notation (13C, 57Fe, Mo-92).  Raises ValueError on invalid input.
    """
    formula_str = formula_str.strip()
    if not formula_str:
        raise ValueError("Formula string is empty.")

    # "Mo-92" or "Fe-56"  →  "[92Mo]" / "[56Fe]"
    dash_match = re.match(r'^([A-Z][a-z]?)-(\d+)$', formula_str)
    if dash_match:
        symbol, mass_num = dash_match.groups()
        formula_str = f"[{mass_num}{symbol}]"

    # "92Mo"  →  "[92Mo]"
    elif re.match(r'^(\d+)([A-Z][a-z]?)$', formula_str):
        specific = re.match(r'^(\d+)([A-Z][a-z]?)$', formula_str)
        mass_num, symbol = specific.groups()
        formula_str = f"[{mass_num}{symbol}]"

    try:
        return Formula(formula_str)
    except Exception as exc:
        raise ValueError(f"Could not parse formula '{formula_str}': {exc}") from exc


def isotope_rows(
    formula_str: str,
    double_charge: bool = False,
    min_fraction: float = 1e-6,
) -> pd.DataFrame:
    """
    Generate isotope table rows for a given formula.

    Returns a DataFrame with columns:
        label, formula, nominal_mass, mz, abundance,
        relative_intensity, charge
    """
    f = _parse_formula(formula_str)
    spectrum = f.spectrum()

    # For specific single-isotope inputs (e.g. "Mo-95", "95Mo"), molmass
    # normalises fraction=1.0.  Look up real values from the parent element.
    parent_spectrum = None
    _is_specific_isotope = bool(
        re.match(r'^\d+[A-Z]', formula_str.strip()) or
        re.match(r'^[A-Z][a-z]?-\d+$', formula_str.strip()) or
        re.match(r'^\[[^\]]+\]$', formula_str.strip())
    )
    if len(spectrum) == 1 and _is_specific_isotope:
        clean_sym = re.sub(r'[\[\]\d]', '', formula_str.strip())
        clean_sym = re.sub(r'^([A-Z][a-z]?).*$', r'\1', clean_sym)
        try:
            parent_spectrum = Formula(clean_sym).spectrum()
        except Exception:
            parent_spectrum = None

    rows = []
    for iso in spectrum.values():
        if iso.fraction < min_fraction:
            continue

        charge = 2 if double_charge else 1
        mz     = iso.mz / charge

        if parent_spectrum is not None and iso.massnumber in parent_spectrum:
            real_iso           = parent_spectrum[iso.massnumber]
            abundance          = float(real_iso.fraction)
            relative_intensity = float(real_iso.intensity)
        else:
            abundance          = float(iso.fraction)
            relative_intensity = float(iso.intensity)

        # Build human-readable label
        clean = re.sub(r'[\[\]]', '', formula_str.strip())
        clean = re.sub(r'^([A-Z][a-z]?)-\d+$', r'\1', clean)
        clean = re.sub(r'^\d+([A-Z][a-z]?)$',  r'\1', clean)
        label = f"{clean}-{iso.massnumber}" + (" (2+)" if double_charge else "")

        rows.append({
            "label":              label,
            "formula":            formula_str.strip(),
            "nominal_mass":       int(iso.massnumber) / charge,
            "mz":                 float(mz),
            "abundance":          abundance,
            "relative_intensity": relative_intensity,
            "charge":             charge,
        })

    if not rows:
        raise ValueError(
            f"No isotopes found for '{formula_str}' above min_fraction={min_fraction}."
        )

    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 2 — H5 lightweight reader
# ═════════════════════════════════════════════════════════════════════════════

def _safe_float(val) -> float | None:
    """Return float(val) or None if val is None / NaN / uncastable."""
    try:
        v = float(val)
        return None if np.isnan(v) else v
    except (TypeError, ValueError):
        return None


def _read_h5_light(file_bytes: bytes) -> dict:
    """
    Read only SumSpectrum + calibration coefficients from an H5 file.
    Never touches TofData, so even huge files load quickly.

    Reads index-domain params for calibration use, plus time-domain
    params for display only.

    Returns dict with keys:
        sum_spectrum, record_size,
        p1_orig, p2_orig,             — mass_calibration_p1/p2 (time-domain, for display)
        p1_recalc, p2_recalc          — recalculated index-domain params (for calibration)
        p1_orig_time, p2_orig_time    — time-domain original (for display)
        p1_recalc_time, p2_recalc_time — time-domain recalculated (for display)
    """
    result = dict(
        sum_spectrum=None,
        record_size=None,
        p1_orig=None,
        p2_orig=None,
        p1_recalc=None,
        p2_recalc=None,
        p1_orig_time=None,
        p2_orig_time=None,
        p1_recalc_time=None,
        p2_recalc_time=None,
    )
    with h5py.File(io.BytesIO(file_bytes), "r") as f:
        if "FullSpectra" in f and "SumSpectrum" in f["FullSpectra"]:
            result["sum_spectrum"] = f["FullSpectra"]["SumSpectrum"][:].astype(np.float64)
        if "Descriptor" in f:
            result["record_size"] = int(f["Descriptor"].attrs.get("recordSize", 0))
        if "Metadata" in f:
            attrs = dict(f["Metadata"].attrs)
            # Time-domain originals (what the file stores as "original")
            result["p1_orig_time"] = _safe_float(attrs.get("mass_calibration_p1"))
            result["p2_orig_time"] = _safe_float(attrs.get("mass_calibration_p2"))
            # Index-domain recalculated (used for actual calibration)
            result["p1_recalc"]    = _safe_float(attrs.get("mass_calibration_recalculated_p1"))
            result["p2_recalc"]    = _safe_float(attrs.get("mass_calibration_recalculated_p2"))
            # Time-domain recalculated (for display)
            result["p1_recalc_time"] = _safe_float(attrs.get("mass_calibration_recalculated_p1_time"))
            result["p2_recalc_time"] = _safe_float(attrs.get("mass_calibration_recalculated_p2_time"))
            # p1_orig/p2_orig: prefer recalculated index, fall back to converting time-domain
            result["p1_orig"] = result["p1_recalc"]
            result["p2_orig"] = result["p2_recalc"]
    return result


def _apply_calibration(indices: np.ndarray, p1: float, p2: float) -> np.ndarray:
    """
    Convert sample indices to m/z using the TOF calibration model:
        index = p1 * sqrt(mass) + p2
    Inverted:
        mass = ((index - p2) / p1) ** 2
    p1 and p2 here are the index-domain parameters stored in /Metadata.
    """
    with np.errstate(invalid="ignore"):
        mz = ((indices - p2) / p1) ** 2
    return np.where(np.isfinite(mz) & (mz >= 0), mz, 0.0)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Gaussian plotting helpers
# ═════════════════════════════════════════════════════════════════════════════

_GAUSSIAN_POINTS = 400

_COLOURS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]


def _hex_to_rgba(hex_color: str, alpha: float) -> str:
    hex_color = hex_color.lstrip("#")
    r, g, b = (int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    return f"rgba({r},{g},{b},{alpha})"


def _gaussian_y(x: np.ndarray, amp: float, centre: float, sigma: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - centre) / sigma) ** 2)


def _make_gaussian_trace(
    centre: float,
    amplitude: float,
    sigma: float,
    colour: str,
    name: str,
    showlegend: bool,
    legendgroup: str,
    hover_text: str,
) -> go.Scatter:
    half_width = 4.5 * sigma
    x = np.linspace(centre - half_width, centre + half_width, _GAUSSIAN_POINTS)
    y = _gaussian_y(x, amplitude, centre, sigma)

    return go.Scatter(
        x=x,
        y=y,
        mode="lines",
        fill="tozeroy",
        fillcolor=_hex_to_rgba(colour, 0.20),
        line=dict(color=colour, width=1.8),
        name=name,
        legendgroup=legendgroup,
        showlegend=showlegend,
        hovertemplate=(
            f"<b>{name}</b><br>"
            "m/z: %{x:.4f}<br>"
            "Intensity: %{y:.4f}<br>"
            f"{hover_text}"
            "<extra></extra>"
        ),
    )


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 4 — Session-state initialisation
# ═════════════════════════════════════════════════════════════════════════════

def _init_state():
    defaults = {
        "active_tool":            "Isotope Plotter",  # which tool is shown in the main panel
        "scaling_factors":        {},       # formula → float
        "dc_scaling_factors":     {},       # formula → float (for 2+ peaks)
        "double_charge_families": set(),    # set of formula strings in 2+ mode
        "peak_width":             0.025,    # Gaussian σ in Da
        "p1":                     None,
        "p2":                     None,
        "spectrum_data":          None,     # (mz_array, intensity_array)
        "loaded_filename":        None,     # name of currently loaded file
        "formula_list":           [],  # persistent ordered list
        "cal_points":             pd.DataFrame({"Peak Index": [0.0, 0.0], "Known m/z": ["", ""]}),
        "cal_sampling_rate":      2.0,
        "cal_fitted_params":      None,   # (p1_idx, p2_idx, p1_time, p2_time)
        "cal_original_params":    None,   # (p1_idx, p2_idx) from file — never overwritten
        "cal_orig_time":          None,   # (p1_time, p2_time) from file — for display only
        "raw_spectrum":           None,   # (indices, raw_intensity) — never baseline-subtracted
    }
    for k, v in defaults.items():
        if k not in st.session_state:
            st.session_state[k] = v


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 5 — Build combined isotope DataFrame
# ═════════════════════════════════════════════════════════════════════════════

def _build_iso_df(formulas: list[str]) -> pd.DataFrame:
    rows = []
    for formula in formulas:
        dc = formula in st.session_state.double_charge_families
        try:
            df_1 = isotope_rows(formula, double_charge=False)
        except ValueError as e:
            st.warning(f"⚠️ {e}")
            continue

        sf     = st.session_state.scaling_factors.get(formula, 0.5)
        sf_dc  = st.session_state.dc_scaling_factors.get(formula, 0.5)
        max_ri = df_1["relative_intensity"].max()

        for _, r in df_1.iterrows():
            scaled    = (r["relative_intensity"] / max_ri) * sf    if max_ri > 0 else sf
            scaled_dc = (r["relative_intensity"] / max_ri) * sf_dc if max_ri > 0 else sf_dc
            rows.append({
                "Formula":        formula,
                "Label":          r["label"],
                "m/z":            round(r["mz"], 5),
                "m/z (2+)":       round(r["mz"] / 2, 5),
                "Abundance (%)":  round(r["abundance"] * 100, 3),
                "Rel. Intensity": round(r["relative_intensity"], 4),
                "2+":             dc,
                "_scaled":        scaled,
                "_scaled_dc":     scaled_dc,
                "_charge":        1,
            })
    return pd.DataFrame(rows)


# ═════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Streamlit app
# ═════════════════════════════════════════════════════════════════════════════

st.set_page_config(page_title="Exum – Mass Spectra Tools", layout="wide")
st.title("Mass Spectra Tools")

_init_state()

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Mass Spectra Tools")
    st.session_state.active_tool = st.radio(
        "Select a tool",
        ["Isotope Plotter", "m/z Lookup", "Oxide Converter"],
        index=["Isotope Plotter", "m/z Lookup", "Oxide Converter"].index(st.session_state.active_tool),
        key="active_tool_radio",
    )
    st.divider()

    if st.session_state.active_tool == "Isotope Plotter":
        st.header("Formulas")
        add_input = st.text_input(
            "Add formula(s)",
            "",
            placeholder="e.g. Mo, Fe2O3, TiO, Ni",
            help="Type element or molecular formulas (e.g. Fe, TiO, Fe2O3). Specific isotopes (57Fe, Mo-92) are not supported — enter the element family instead.",
            key="formula_add_input",
        )

        if st.button("Add", type="primary", key="formula_add_btn", width="stretch"):
            new_formulas = [f.strip() for f in add_input.split(",") if f.strip()]
            added = False
            rejected = []
            for nf in new_formulas:
                if _is_single_isotope(nf):
                    rejected.append(nf)
                elif nf not in st.session_state.formula_list:
                    st.session_state.formula_list.append(nf)
                    added = True
            if rejected:
                st.warning(
                    f"**{', '.join(rejected)}** — specific isotopes are not allowed. "
                    f"Enter an element or molecular formula (e.g. `Fe`, `TiO`, `Fe2O3`) "
                    f"to plot the full family."
                )
            if not new_formulas:
                st.warning("Enter at least one formula first.")
            elif not added and not rejected:
                st.info("Formula(s) already in the list.")
            # No st.rerun() — button click already triggers a full rerun

    if st.session_state.active_tool == "m/z Lookup":
        #NEW CODE 1
        st.subheader("Reverse Lookup (m/z → isotopes)")

        target_mz = st.number_input(
            "Enter m/z",
            min_value=0.0,
            value=0.0,
            step=0.1,
            key="reverse_mz_input"
        )
        element_filter = st.text_input(
            "Filter by element(s)",
            "",
            placeholder="e.g. Fe, Ti, Si, Li",
            help="Leave blank to search all species. Enter one or more element symbols separated by commas.",
        )

        tolerance = st.number_input(
            "Tolerance (Da)",
            min_value=0.0001,
            value=0.1,
            step=0.01,
            key="reverse_tol_input"
        )

        run_reverse = st.button("Find Matches")

        #END NEW CODE 1
    
    if st.session_state.active_tool == "Oxide Converter":

        st.header("Oxide → Element Converter")

        uploaded_file = st.file_uploader(
            "Upload oxide standard CSV",
            type=["csv"],
            key="oxide_converter_upload",
        )
            
    
    formulas = st.session_state.formula_list

    if st.session_state.active_tool == "Isotope Plotter":
        st.divider()
        st.header("Peak Appearance")
        peak_width = st.slider(
            "Peak width σ (Da)", min_value=0.005, max_value=0.3,
            value=st.session_state.peak_width, step=0.005, format="%.3f",
        )
        st.session_state.peak_width = peak_width

        st.divider()
        st.header("Upload Spectrum")
        uploaded_file = st.file_uploader(
            "Mass spectrum (.h5 or .csv)", type=["h5", "csv"],
            help="H5: only SumSpectrum + calibration coefficients are read (memory-safe).",
        )

        if uploaded_file is None:
            # Only clear if we previously had a file loaded — not on every rerun
            if st.session_state.loaded_filename is not None:
                st.session_state.spectrum_data       = None
                st.session_state.raw_spectrum        = None
                st.session_state.p1                  = None
                st.session_state.p2                  = None
                st.session_state.cal_original_params = None
                st.session_state.cal_orig_time       = None
                st.session_state.cal_fitted_params   = None
                st.session_state.cal_using_new       = False
                st.session_state._cal_residuals      = None
                st.session_state.loaded_filename     = None
                st.rerun()
        else:
            file_bytes = uploaded_file.read()
            fname      = uploaded_file.name

            if fname.endswith(".csv"):
                try:
                    csv_df  = pd.read_csv(io.BytesIO(file_bytes), header=0)
                    mz_arr  = csv_df.iloc[:, 0].values.astype(np.float64)
                    int_arr = csv_df.iloc[:, 1].values.astype(np.float64)
                    st.session_state.spectrum_data = (mz_arr, int_arr)
                    st.session_state.p1 = None
                    st.session_state.p2 = None
                    st.session_state.loaded_filename = fname
                    st.success("CSV loaded.")
                except Exception as e:
                    st.error(f"CSV read error: {e}")

            elif fname.endswith(".h5"):
                if not _H5_AVAILABLE:
                    st.error("h5py is not installed.  Run: pip install h5py")
                else:
                    try:
                        h5   = _read_h5_light(file_bytes)
                        spec = h5["sum_spectrum"]
                        if spec is None:
                            st.error("No SumSpectrum found in file.")
                        else:
                            n       = h5["record_size"] or len(spec)
                            indices = np.arange(n, dtype=np.float64)

                            # Prefer recalculated coefficients; fall back to originals
                            p1_use = h5["p1_recalc"] if h5["p1_recalc"] is not None else h5["p1_orig"]
                            p2_use = h5["p2_recalc"] if h5["p2_recalc"] is not None else h5["p2_orig"]

                            if p1_use is not None and p2_use is not None:
                                st.session_state.p1 = p1_use
                                st.session_state.p2 = p2_use
                                mz_arr = _apply_calibration(indices, p1_use, p2_use)
                            else:
                                st.warning("No calibration coefficients found — using sample indices as x-axis.")
                                mz_arr = indices
                                st.session_state.p1 = None
                                st.session_state.p2 = None

                            # Baseline correction (optional)
                            if _BASELINE_AVAILABLE:
                                bl      = Baseline(x_data=mz_arr).snip(spec, max_half_window=18, decreasing=True)[0]
                                int_arr = np.clip(spec - bl, 0, None)
                            else:
                                int_arr = np.clip(spec, 0, None)

                            st.session_state.spectrum_data = (mz_arr, int_arr)
                            # Also keep raw index-domain data for the calibration tool
                            st.session_state.raw_spectrum = (indices, spec)
                            # Store original calibration params once for revert
                            if p1_use is not None and p2_use is not None:
                                st.session_state.cal_original_params = (p1_use, p2_use)
                            # Store time-domain originals for display only
                            st.session_state.cal_orig_time = (
                                h5["p1_orig_time"], h5["p2_orig_time"]
                            ) if h5["p1_orig_time"] is not None else None
                            st.session_state.loaded_filename = fname
                            st.success(f"H5 loaded — {len(spec):,} points.")
                    except Exception as e:
                        st.error(f"H5 read error: {e}")

        # (calibration UI moved to main panel)

if st.session_state.active_tool == "Isotope Plotter":
    # ── Main panel ────────────────────────────────────────────────────────────────

    if not formulas:
        st.info("Add chemical formulas in the sidebar to begin.")
        # Still show the spectrum if a file is loaded
        if st.session_state.spectrum_data is not None:
            mz_arr, int_arr = st.session_state.spectrum_data
            fig_bare = go.Figure()
            fig_bare.add_trace(go.Scatter(
                x=mz_arr, y=int_arr,
                mode="lines",
                name="Mass Spectrum",
                line=dict(color="black", width=1.2),
            ))
            fig_bare.update_layout(
                xaxis_title="m/z",
                yaxis_title="Intensity",
                autosize=True,
                margin=dict(t=40, b=60, l=60, r=20),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
            )
            fig_bare.update_xaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
            fig_bare.update_yaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
            st.plotly_chart(fig_bare)

        # ── Mass Calibration expander ─────────────────────────────────────────────────
        if st.session_state.raw_spectrum is not None:
            with st.expander("⚖️ Mass Calibration", expanded=False):
                from scipy.optimize import curve_fit as _curve_fit
                import re as _re

                # ── TOF model ────────────────────────────────────────────────────────
                def _tof_model(mass, p1, p2):
                    """index = p1 * sqrt(mass) + p2"""
                    return p1 * np.sqrt(mass) + p2

                # ── Formula → exact m/z resolver (copied from 1_Processing.py) ──────
                def _resolve_label_to_mz(raw: str) -> tuple[float, str]:
                    """
                    Resolve a label string to an exact m/z value.

                    Accepts:
                      - Plain numbers:           "63.9291"
                      - Single isotopes:         "Cu-63", "63Cu", "[63Cu]"
                      - Molecular formulas:      "TiO", "Si2"  (only if single isotope result)
                      - Molecular + mass number: "TiO-64", "Si2-56"
                        → looks up the peak at that nominal mass in the molecule spectrum
                      - Leading-mass molecules:  "64TiO"
                        → converted to "TiO-64" first
                    Returns (mz_float, "") on success or (nan, error_str) on failure.
                    """
                    s = raw.strip()

                    # Plain number
                    try:
                        return float(s), ""
                    except ValueError:
                        pass

                    # Leading-mass molecular notation: "64TiO" → "TiO-64"
                    _lm = _re.match(r'^(\d+)([A-Z][A-Za-z0-9]+)$', s)
                    if _lm:
                        s = f"{_lm.group(2)}-{_lm.group(1)}"

                    # "Formula-MassNumber" notation: "TiO-64", "Cu-63"
                    _fm = _re.match(r'^([A-Za-z][A-Za-z0-9]*)-(\d+)$', s)
                    if _fm:
                        _formula_part = _fm.group(1)
                        _target_mass  = int(_fm.group(2))
                        try:
                            _f    = _parse_formula(_formula_part)
                            _spec = _f.spectrum()
                            if len(_spec) == 1:
                                return float(next(iter(_spec.values())).mz), ""
                            if _target_mass in _spec:
                                return float(_spec[_target_mass].mz), ""
                            _avail = sorted(_spec.keys())
                            return float("nan"), (
                                f"Mass {_target_mass} not found in {_formula_part} spectrum. "
                                f"Available: {', '.join(str(m) for m in _avail)}"
                            )
                        except Exception as _e:
                            return float("nan"), str(_e)

                    # Plain formula (no mass number)
                    try:
                        _f    = _parse_formula(s)
                        _spec = _f.spectrum()
                        if len(_spec) == 1:
                            return float(next(iter(_spec.values())).mz), ""
                        _base         = s.translate(str.maketrans("", "", "[]0123456789")).strip("-")
                        _example_mass = max(_spec.values(), key=lambda _i: _i.fraction).massnumber
                        return float("nan"), (
                            f"**{raw}** matches multiple isotopes — "
                            f"be specific, e.g. `{_base}-{_example_mass}`."
                        )
                    except Exception as _e:
                        return float("nan"), str(_e)

                st.caption(
                    "Enter an **Index Position** (hover the spectrum below to read it) and a "
                    "**Known m/z** label. Labels can be plain numbers (`63.9291`), isotopes "
                    "(`Cu-63`, `63Cu`), or molecular formulas with a mass number "
                    "(`TiO-64`, `64TiO`, `Si2-56`). Press **Run Calibration** to fit."
                )

                # ── Calibration points table ──────────────────────────────────────────
                st.markdown("**Calibration Points**")
                # Key design (matching 1_Processing.py): seed the editor ONCE from a
                # fixed default DataFrame and never overwrite it from the return value.
                # Streamlit keeps the widget's own edited state under the key — writing
                # back on every rerun is what causes the double-click reset.
                _CAL_SEED = pd.DataFrame({
                    "Peak Index": pd.array([None, None], dtype="object"),
                    "Known m/z":  ["", ""],
                })
                cal_points_raw = st.data_editor(
                    _CAL_SEED,
                    num_rows="dynamic",
                    width="stretch",
                    column_config={
                        "Peak Index": st.column_config.NumberColumn(
                            "Index Position", format="%.1f",
                        ),
                        "Known m/z": st.column_config.TextColumn(
                            "Known m/z",
                            help=(
                                "Plain number (e.g. 63.9291), isotope (Cu-63, 63Cu), "
                                "or formula+mass (TiO-64, 64TiO, Si2-56)."
                            ),
                        ),
                    },
                    key="cal_table_editor",
                )

                # ── Resolve labels → numeric m/z live ────────────────────────────────
                _resolved_mz    = []
                _resolve_errors = []
                for _ri, _row in cal_points_raw.iterrows():
                    _raw = str(_row.get("Known m/z", "") or "").strip()
                    if not _raw:
                        _resolved_mz.append(float("nan"))
                        continue
                    _mzv, _err = _resolve_label_to_mz(_raw)
                    _resolved_mz.append(_mzv)
                    if _err:
                        _resolve_errors.append(f"Row {_ri + 1}: {_err}")

                for _err in _resolve_errors:
                    st.warning(_err)

                # Live preview of resolved labels
                _preview_rows = []
                for _ri, _row in cal_points_raw.iterrows():
                    _raw = str(_row.get("Known m/z", "") or "").strip()
                    _mzv = _resolved_mz[_ri] if _ri < len(_resolved_mz) else float("nan")
                    if _raw and pd.notna(_mzv):
                        try:
                            float(_raw)   # skip plain numbers — no need to preview
                        except ValueError:
                            _preview_rows.append(f"`{_raw}` → **{_mzv:.6f}**")
                if _preview_rows:
                    st.caption("Resolved: " + " · ".join(_preview_rows))

                # Build numeric cal_table for fitting
                cal_table = pd.DataFrame({
                    "Peak Index": cal_points_raw["Peak Index"].values,
                    "Known m/z":  _resolved_mz,
                })
                cal_valid = cal_table.dropna(subset=["Peak Index", "Known m/z"])
                cal_valid = cal_valid[(cal_valid["Peak Index"] > 0) & (cal_valid["Known m/z"] > 0)]
                cal_ready = len(cal_valid) >= 2 and not _resolve_errors
                if not cal_ready and len(cal_valid) < 2:
                    st.caption("Enter at least 2 valid index / m/z pairs to run calibration.")

                sampling_rate = st.number_input(
                    "Digitiser sampling rate (GS/s)",
                    min_value=0.01, max_value=100.0,
                    value=st.session_state.cal_sampling_rate, step=0.1,
                    help="Converts index-domain p-values to time-domain (s) p-values.",
                    key="cal_sampling_rate_input",
                )
                st.session_state.cal_sampling_rate = sampling_rate

                # ── Buttons ───────────────────────────────────────────────────────────
                _has_orig  = st.session_state.cal_original_params is not None
                _has_new   = st.session_state.cal_fitted_params is not None
                # Use an explicit flag set when we apply/revert — avoids float equality fragility
                if "cal_using_new" not in st.session_state:
                    st.session_state.cal_using_new = False
                _using_new = st.session_state.cal_using_new

                col_run, col_revert = st.columns([1, 1])
                with col_run:
                    run_cal = st.button(
                        "Run & Apply Calibration", type="primary", key="cal_run_btn",
                        disabled=not cal_ready,
                    )
                with col_revert:
                    revert_cal = st.button(
                        "↩ Revert to Original", key="cal_revert_btn",
                        disabled=not (_has_orig and _has_new),
                        help="Switch back to the calibration parameters from the file.",
                    )

                if run_cal and cal_ready:
                    masses_fit  = cal_valid["Known m/z"].values
                    indices_fit = cal_valid["Peak Index"].values
                    try:
                        import warnings as _warnings
                        from scipy.optimize import OptimizeWarning as _OptimizeWarning
                        with _warnings.catch_warnings():
                            _warnings.simplefilter("ignore", _OptimizeWarning)
                            popt, _    = _curve_fit(_tof_model, masses_fit, indices_fit, p0=[100, 0])
                            p1_i, p2_i = popt

                            times_s    = indices_fit / (sampling_rate * 1e9)
                            popt_t, _  = _curve_fit(_tof_model, masses_fit, times_s, p0=[1e-7, 0])
                            p1_t, p2_t = popt_t

                        st.session_state.cal_fitted_params = (p1_i, p2_i, p1_t, p2_t)

                        # Residuals
                        fitted_idx = _tof_model(masses_fit, p1_i, p2_i)
                        st.session_state._cal_residuals = [
                            (float(m), float(ii), float(fi), float(ii - fi))
                            for m, ii, fi in zip(masses_fit, indices_fit, fitted_idx)
                        ]

                        # Apply immediately
                        _, raw_int = st.session_state.raw_spectrum
                        n_pts  = len(raw_int)
                        new_mz = _apply_calibration(np.arange(n_pts, dtype=np.float64), p1_i, p2_i)
                        if _BASELINE_AVAILABLE:
                            bl      = Baseline(x_data=new_mz).snip(raw_int, max_half_window=18, decreasing=True)[0]
                            new_int = np.clip(raw_int - bl, 0, None)
                        else:
                            new_int = np.clip(raw_int, 0, None)
                        st.session_state.spectrum_data = (new_mz, new_int)
                        st.session_state.p1 = p1_i
                        st.session_state.p2 = p2_i
                        st.success("Calibration fitted and applied to spectrum.")

                    except Exception as _e:
                        st.error(f"Curve fitting failed: {_e}")

                if revert_cal and st.session_state.cal_original_params is not None:
                    p1_o, p2_o = st.session_state.cal_original_params
                    _, raw_int = st.session_state.raw_spectrum
                    n_pts  = len(raw_int)
                    rev_mz = _apply_calibration(np.arange(n_pts, dtype=np.float64), p1_o, p2_o)
                    if _BASELINE_AVAILABLE:
                        bl      = Baseline(x_data=rev_mz).snip(raw_int, max_half_window=18, decreasing=True)[0]
                        rev_int = np.clip(raw_int - bl, 0, None)
                    else:
                        rev_int = np.clip(raw_int, 0, None)
                    st.session_state.spectrum_data = (rev_mz, rev_int)
                    st.session_state.p1 = p1_o
                    st.session_state.p2 = p2_o
                    st.session_state.cal_using_new = False
                    st.rerun()

                # ── Display parameters (original vs fitted) ──────────────────────────
                if st.session_state.cal_fitted_params is not None or st.session_state.cal_original_params is not None:
                    st.markdown("**Calibration Parameters** *(not written to file)*")
                    st.caption(f"Currently active: {'**New (fitted)**' if _using_new else '**Original (from file)**'}")
                    _pc1, _pc2, _pc3, _pc4 = st.columns(4)
                    if st.session_state.cal_orig_time is not None:
                        p1_ot, p2_ot = st.session_state.cal_orig_time
                        with _pc1:
                            st.markdown("*Original — Time (from file)*")
                            st.metric("p1 (s)", f"{p1_ot:.6e}")
                            st.metric("p2 (s)", f"{p2_ot:.6e}")
                    if st.session_state.cal_fitted_params is not None:
                        p1_i, p2_i, p1_t, p2_t = st.session_state.cal_fitted_params
                        with _pc2:
                            st.markdown("*New — Index*")
                            st.metric("p1", f"{p1_i:.6e}")
                            st.metric("p2", f"{p2_i:.6e}")
                        with _pc3:
                            st.markdown("*New — Time*")
                            _d1t = p1_t - p1_ot if st.session_state.cal_orig_time else None
                            _d2t = p2_t - p2_ot if st.session_state.cal_orig_time else None
                            st.metric("p1 (s)", f"{p1_t:.6e}", delta=f"{_d1t:+.3e}" if _d1t is not None else None)
                            st.metric("p2 (s)", f"{p2_t:.6e}", delta=f"{_d2t:+.3e}" if _d2t is not None else None)

                    if st.session_state.get("_cal_residuals"):
                        st.markdown("**Residuals**")
                        res_df = pd.DataFrame(
                            st.session_state._cal_residuals,
                            columns=["Known m/z", "Measured Index", "Fitted Index", "Residual (samples)"],
                        )
                        st.dataframe(
                            res_df.style.format({
                                "Known m/z":           "{:.6f}",
                                "Measured Index":       "{:.1f}",
                                "Fitted Index":         "{:.2f}",
                                "Residual (samples)":   "{:+.3f}",
                            }),
                            hide_index=True,
                            width="stretch",
                        )

                # ── Index-domain spectrum viewer ──────────────────────────────────────
                st.markdown("**Raw Spectrum (index domain)** — hover to read peak indices")
                idx_arr, raw_int = st.session_state.raw_spectrum
                fig_cal = go.Figure()
                fig_cal.add_trace(go.Scatter(
                    x=idx_arr, y=raw_int,
                    mode="lines",
                    line=dict(color="#1f77b4", width=1),
                    name="Raw",
                    hovertemplate="Index: %{x:.0f}<br>Intensity: %{y:.2f}<extra></extra>",
                ))
                fig_cal.update_layout(
                    xaxis_title="Sample Index",
                    yaxis_title="Intensity",
                    autosize=True,
                    margin=dict(t=20, b=40, l=60, r=20),
                    plot_bgcolor="rgba(0,0,0,0)",
                    paper_bgcolor="rgba(0,0,0,0)",
                    showlegend=False,
                )
                fig_cal.update_xaxes(showgrid=False, zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
                fig_cal.update_yaxes(showgrid=False, zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
                st.plotly_chart(fig_cal)

        st.stop()


    # ── Mass Calibration expander ─────────────────────────────────────────────────
    if st.session_state.raw_spectrum is not None:
        with st.expander("⚖️ Mass Calibration", expanded=False):
            from scipy.optimize import curve_fit as _curve_fit
            import re as _re

            # ── TOF model ────────────────────────────────────────────────────────
            def _tof_model(mass, p1, p2):
                """index = p1 * sqrt(mass) + p2"""
                return p1 * np.sqrt(mass) + p2

            # ── Formula → exact m/z resolver (copied from 1_Processing.py) ──────
            def _resolve_label_to_mz(raw: str) -> tuple[float, str]:
                """
                Resolve a label string to an exact m/z value.

                Accepts:
                  - Plain numbers:           "63.9291"
                  - Single isotopes:         "Cu-63", "63Cu", "[63Cu]"
                  - Molecular formulas:      "TiO", "Si2"  (only if single isotope result)
                  - Molecular + mass number: "TiO-64", "Si2-56"
                    → looks up the peak at that nominal mass in the molecule spectrum
                  - Leading-mass molecules:  "64TiO"
                    → converted to "TiO-64" first
                Returns (mz_float, "") on success or (nan, error_str) on failure.
                """
                s = raw.strip()

                # Plain number
                try:
                    return float(s), ""
                except ValueError:
                    pass

                # Leading-mass molecular notation: "64TiO" → "TiO-64"
                _lm = _re.match(r'^(\d+)([A-Z][A-Za-z0-9]+)$', s)
                if _lm:
                    s = f"{_lm.group(2)}-{_lm.group(1)}"

                # "Formula-MassNumber" notation: "TiO-64", "Cu-63"
                _fm = _re.match(r'^([A-Za-z][A-Za-z0-9]*)-(\d+)$', s)
                if _fm:
                    _formula_part = _fm.group(1)
                    _target_mass  = int(_fm.group(2))
                    try:
                        _f    = _parse_formula(_formula_part)
                        _spec = _f.spectrum()
                        if len(_spec) == 1:
                            return float(next(iter(_spec.values())).mz), ""
                        if _target_mass in _spec:
                            return float(_spec[_target_mass].mz), ""
                        _avail = sorted(_spec.keys())
                        return float("nan"), (
                            f"Mass {_target_mass} not found in {_formula_part} spectrum. "
                            f"Available: {', '.join(str(m) for m in _avail)}"
                        )
                    except Exception as _e:
                        return float("nan"), str(_e)

                # Plain formula (no mass number)
                try:
                    _f    = _parse_formula(s)
                    _spec = _f.spectrum()
                    if len(_spec) == 1:
                        return float(next(iter(_spec.values())).mz), ""
                    _base         = s.translate(str.maketrans("", "", "[]0123456789")).strip("-")
                    _example_mass = max(_spec.values(), key=lambda _i: _i.fraction).massnumber
                    return float("nan"), (
                        f"**{raw}** matches multiple isotopes — "
                        f"be specific, e.g. `{_base}-{_example_mass}`."
                    )
                except Exception as _e:
                    return float("nan"), str(_e)

            st.caption(
                "Enter an **Index Position** (hover the spectrum below to read it) and a "
                "**Known m/z** label. Labels can be plain numbers (`63.9291`), isotopes "
                "(`Cu-63`, `63Cu`), or molecular formulas with a mass number "
                "(`TiO-64`, `64TiO`, `Si2-56`). Press **Run Calibration** to fit."
            )

            # ── Calibration points table ──────────────────────────────────────────
            st.markdown("**Calibration Points**")
            # Key design (matching 1_Processing.py): seed the editor ONCE from a
            # fixed default DataFrame and never overwrite it from the return value.
            # Streamlit keeps the widget's own edited state under the key — writing
            # back on every rerun is what causes the double-click reset.
            _CAL_SEED = pd.DataFrame({
                "Peak Index": pd.array([None, None], dtype="object"),
                "Known m/z":  ["", ""],
            })
            cal_points_raw = st.data_editor(
                _CAL_SEED,
                num_rows="dynamic",
                width="stretch",
                column_config={
                    "Peak Index": st.column_config.NumberColumn(
                        "Index Position", format="%.1f",
                    ),
                    "Known m/z": st.column_config.TextColumn(
                        "Known m/z",
                        help=(
                            "Plain number (e.g. 63.9291), isotope (Cu-63, 63Cu), "
                            "or formula+mass (TiO-64, 64TiO, Si2-56)."
                        ),
                    ),
                },
                key="cal_table_editor",
            )

            # ── Resolve labels → numeric m/z live ────────────────────────────────
            _resolved_mz    = []
            _resolve_errors = []
            for _ri, _row in cal_points_raw.iterrows():
                _raw = str(_row.get("Known m/z", "") or "").strip()
                if not _raw:
                    _resolved_mz.append(float("nan"))
                    continue
                _mzv, _err = _resolve_label_to_mz(_raw)
                _resolved_mz.append(_mzv)
                if _err:
                    _resolve_errors.append(f"Row {_ri + 1}: {_err}")

            for _err in _resolve_errors:
                st.warning(_err)

            # Live preview of resolved labels
            _preview_rows = []
            for _ri, _row in cal_points_raw.iterrows():
                _raw = str(_row.get("Known m/z", "") or "").strip()
                _mzv = _resolved_mz[_ri] if _ri < len(_resolved_mz) else float("nan")
                if _raw and pd.notna(_mzv):
                    try:
                        float(_raw)   # skip plain numbers — no need to preview
                    except ValueError:
                        _preview_rows.append(f"`{_raw}` → **{_mzv:.6f}**")
            if _preview_rows:
                st.caption("Resolved: " + " · ".join(_preview_rows))

            # Build numeric cal_table for fitting
            cal_table = pd.DataFrame({
                "Peak Index": cal_points_raw["Peak Index"].values,
                "Known m/z":  _resolved_mz,
            })
            cal_valid = cal_table.dropna(subset=["Peak Index", "Known m/z"])
            cal_valid = cal_valid[(cal_valid["Peak Index"] > 0) & (cal_valid["Known m/z"] > 0)]
            cal_ready = len(cal_valid) >= 2 and not _resolve_errors
            if not cal_ready and len(cal_valid) < 2:
                st.caption("Enter at least 2 valid index / m/z pairs to run calibration.")

            sampling_rate = st.number_input(
                "Digitiser sampling rate (GS/s)",
                min_value=0.01, max_value=100.0,
                value=st.session_state.cal_sampling_rate, step=0.1,
                help="Converts index-domain p-values to time-domain (s) p-values.",
                key="cal_sampling_rate_input",
            )
            st.session_state.cal_sampling_rate = sampling_rate

            # ── Buttons ───────────────────────────────────────────────────────────
            _has_orig  = st.session_state.cal_original_params is not None
            _has_new   = st.session_state.cal_fitted_params is not None
            # Use an explicit flag set when we apply/revert — avoids float equality fragility
            if "cal_using_new" not in st.session_state:
                st.session_state.cal_using_new = False
            _using_new = st.session_state.cal_using_new

            col_run, col_revert = st.columns([1, 1])
            with col_run:
                run_cal = st.button(
                    "Run & Apply Calibration", type="primary", key="cal_run_btn",
                    disabled=not cal_ready,
                )
            with col_revert:
                revert_cal = st.button(
                    "↩ Revert to Original", key="cal_revert_btn",
                    disabled=not (_has_orig and _has_new),
                    help="Switch back to the calibration parameters from the file.",
                )

            if run_cal and cal_ready:
                masses_fit  = cal_valid["Known m/z"].values
                indices_fit = cal_valid["Peak Index"].values
                try:
                    import warnings as _warnings
                    from scipy.optimize import OptimizeWarning as _OptimizeWarning
                    with _warnings.catch_warnings():
                        _warnings.simplefilter("ignore", _OptimizeWarning)
                        popt, _    = _curve_fit(_tof_model, masses_fit, indices_fit, p0=[100, 0])
                        p1_i, p2_i = popt

                        times_s    = indices_fit / (sampling_rate * 1e9)
                        popt_t, _  = _curve_fit(_tof_model, masses_fit, times_s, p0=[1e-7, 0])
                        p1_t, p2_t = popt_t

                    st.session_state.cal_fitted_params = (p1_i, p2_i, p1_t, p2_t)

                    # Residuals
                    fitted_idx = _tof_model(masses_fit, p1_i, p2_i)
                    st.session_state._cal_residuals = [
                        (float(m), float(ii), float(fi), float(ii - fi))
                        for m, ii, fi in zip(masses_fit, indices_fit, fitted_idx)
                    ]

                    # Apply immediately
                    _, raw_int = st.session_state.raw_spectrum
                    n_pts  = len(raw_int)
                    new_mz = _apply_calibration(np.arange(n_pts, dtype=np.float64), p1_i, p2_i)
                    if _BASELINE_AVAILABLE:
                        bl      = Baseline(x_data=new_mz).snip(raw_int, max_half_window=18, decreasing=True)[0]
                        new_int = np.clip(raw_int - bl, 0, None)
                    else:
                        new_int = np.clip(raw_int, 0, None)
                    st.session_state.spectrum_data = (new_mz, new_int)
                    st.session_state.p1 = p1_i
                    st.session_state.p2 = p2_i
                    st.success("Calibration fitted and applied to spectrum.")

                except Exception as _e:
                    st.error(f"Curve fitting failed: {_e}")

            if revert_cal and st.session_state.cal_original_params is not None:
                p1_o, p2_o = st.session_state.cal_original_params
                _, raw_int = st.session_state.raw_spectrum
                n_pts  = len(raw_int)
                rev_mz = _apply_calibration(np.arange(n_pts, dtype=np.float64), p1_o, p2_o)
                if _BASELINE_AVAILABLE:
                    bl      = Baseline(x_data=rev_mz).snip(raw_int, max_half_window=18, decreasing=True)[0]
                    rev_int = np.clip(raw_int - bl, 0, None)
                else:
                    rev_int = np.clip(raw_int, 0, None)
                st.session_state.spectrum_data = (rev_mz, rev_int)
                st.session_state.p1 = p1_o
                st.session_state.p2 = p2_o
                st.session_state.cal_using_new = False
                st.rerun()

            # ── Display parameters (original vs fitted) ──────────────────────────
            if st.session_state.cal_fitted_params is not None or st.session_state.cal_orig_time is not None:
                st.markdown("**Calibration Parameters** *(not written to file)*")
                st.caption(f"Currently active: {'**New (fitted)**' if _using_new else '**Original (from file)**'}")
                _pc1, _pc2, _pc3, _pc4 = st.columns(4)
                if st.session_state.cal_orig_time is not None:
                    p1_ot, p2_ot = st.session_state.cal_orig_time
                    with _pc1:
                        st.markdown("*Original — Time (from file)*")
                        st.metric("p1 (s)", f"{p1_ot:.6e}")
                        st.metric("p2 (s)", f"{p2_ot:.6e}")
                if st.session_state.cal_fitted_params is not None:
                    p1_i, p2_i, p1_t, p2_t = st.session_state.cal_fitted_params
                    with _pc2:
                        st.markdown("*New — Index*")
                        st.metric("p1", f"{p1_i:.6e}")
                        st.metric("p2", f"{p2_i:.6e}")
                    with _pc3:
                        st.markdown("*New — Time*")
                        _d1t = p1_t - p1_ot if st.session_state.cal_orig_time else None
                        _d2t = p2_t - p2_ot if st.session_state.cal_orig_time else None
                        st.metric("p1 (s)", f"{p1_t:.6e}", delta=f"{_d1t:+.3e}" if _d1t is not None else None)
                        st.metric("p2 (s)", f"{p2_t:.6e}", delta=f"{_d2t:+.3e}" if _d2t is not None else None)

                if st.session_state.get("_cal_residuals"):
                    st.markdown("**Residuals**")
                    res_df = pd.DataFrame(
                        st.session_state._cal_residuals,
                        columns=["Known m/z", "Measured Index", "Fitted Index", "Residual (samples)"],
                    )
                    st.dataframe(
                        res_df.style.format({
                            "Known m/z":           "{:.6f}",
                            "Measured Index":       "{:.1f}",
                            "Fitted Index":         "{:.2f}",
                            "Residual (samples)":   "{:+.3f}",
                        }),
                        hide_index=True,
                        width="stretch",
                    )

            # ── Index-domain spectrum viewer ──────────────────────────────────────
            st.markdown("**Raw Spectrum (index domain)** — hover to read peak indices")
            idx_arr, raw_int = st.session_state.raw_spectrum
            fig_cal = go.Figure()
            fig_cal.add_trace(go.Scatter(
                x=idx_arr, y=raw_int,
                mode="lines",
                line=dict(color="#1f77b4", width=1),
                name="Raw",
                hovertemplate="Index: %{x:.0f}<br>Intensity: %{y:.2f}<extra></extra>",
            ))
            fig_cal.update_layout(
                xaxis_title="Sample Index",
                yaxis_title="Intensity",
                autosize=True,
                margin=dict(t=20, b=40, l=60, r=20),
                plot_bgcolor="rgba(0,0,0,0)",
                paper_bgcolor="rgba(0,0,0,0)",
                showlegend=False,
            )
            fig_cal.update_xaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
            fig_cal.update_yaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
            st.plotly_chart(fig_cal)

    # Ensure default scaling factors exist for any new formulas
    for f in formulas:
        if f not in st.session_state.scaling_factors:
            st.session_state.scaling_factors[f] = 0.5
        if f not in st.session_state.dc_scaling_factors:
            st.session_state.dc_scaling_factors[f] = 0.5

    # Colour map — one colour per formula
    colour_map = {f: _COLOURS[i % len(_COLOURS)] for i, f in enumerate(formulas)}


    # ── Combined isotope controls + detail table ──────────────────────────────────
    st.subheader("Isotope Controls")
    st.caption(
        "Edit **Scale** per row — applies to the whole family.  "
        "Tick **2+** on any isotope in a family to plot both 1+ and 2+ peaks for every isotope in that family."
    )

    # ── Select All / Deselect All ─────────────────────────────────────────────────
    _all_pending = set(formulas) == st.session_state.get("_pending_delete", set())
    _col_sel, _col_desel, _ = st.columns([1, 1, 6])
    with _col_sel:
        if st.button("☑ Select All", key="sel_all_btn", disabled=_all_pending):
            st.session_state._pending_delete = set(formulas)
            st.rerun()
    with _col_desel:
        if st.button("☐ Deselect All", key="desel_all_btn", disabled=not st.session_state.get("_pending_delete")):
            st.session_state._pending_delete = set()
            st.rerun()

    iso_df = _build_iso_df(formulas)

    if iso_df.empty:
        st.error("No valid isotopes could be generated from the entered formulas.")
        st.stop()

    # Stamp current family-level state onto every row before rendering
    iso_df["Scale"]    = iso_df["Formula"].map(st.session_state.scaling_factors)
    iso_df["Scale 2+"] = iso_df["Formula"].map(st.session_state.dc_scaling_factors)
    iso_df["2+"]     = iso_df["Formula"].apply(lambda f: f in st.session_state.double_charge_families)
    iso_df["Delete"] = iso_df["Formula"].apply(lambda f: f in st.session_state.get("_pending_delete", set()))

    edit_cols = ["Formula", "Label", "m/z", "m/z (2+)", "Abundance (%)", "Rel. Intensity", "Scale", "Scale 2+", "2+", "Delete"]

    edited_iso = st.data_editor(
        iso_df[edit_cols],
        column_config={
            "Formula":        st.column_config.TextColumn("Formula",         disabled=True),
            "Label":          st.column_config.TextColumn("Label",           disabled=True),
            "m/z":            st.column_config.NumberColumn("m/z",           disabled=True, format="%.5f"),
            "m/z (2+)":       st.column_config.NumberColumn("m/z (2+)",      disabled=True, format="%.5f"),
            "Abundance (%)":  st.column_config.NumberColumn("Abundance (%)",  disabled=True, format="%.3f"),
            "Rel. Intensity": st.column_config.NumberColumn("Rel. Intensity", disabled=True, format="%.4f"),
            "Scale":    st.column_config.NumberColumn(
                "Scale", min_value=0.0, format="%.6g",
            ),
            "Scale 2+": st.column_config.NumberColumn(
                "Scale 2+", min_value=0.0, format="%.6g",
            ),
            "2+":     st.column_config.CheckboxColumn("2+"),
            "Delete": st.column_config.CheckboxColumn("🗑 Delete"),
        },
        hide_index=True,
        width="stretch",
        key=f"iso_editor_{'_'.join(sorted(formulas))}",
    )

    # Detect changes and propagate family-wide
    # For 2+: if ANY row in a family was toggled, apply that value to the whole family.
    # For scale: any row change updates the whole family.
    changed = False
    for formula in formulas:
        fam_edited = edited_iso[edited_iso["Formula"] == formula]
        fam_orig   = iso_df[iso_df["Formula"] == formula]

        # Scale 1+ — take the first changed value found in the family
        orig_sf = st.session_state.scaling_factors.get(formula, 0.5)
        changed_sf_rows = fam_edited[abs(fam_edited["Scale"] - orig_sf) > 1e-9]
        if not changed_sf_rows.empty:
            st.session_state.scaling_factors[formula] = float(changed_sf_rows.iloc[0]["Scale"])
            changed = True

        # Scale 2+ — separate scale for the double-charge overlay
        orig_sf_dc = st.session_state.dc_scaling_factors.get(formula, 0.5)
        changed_sf_dc_rows = fam_edited[abs(fam_edited["Scale 2+"] - orig_sf_dc) > 1e-9]
        if not changed_sf_dc_rows.empty:
            st.session_state.dc_scaling_factors[formula] = float(changed_sf_dc_rows.iloc[0]["Scale 2+"])
            changed = True

        # 2+ — if any row differs from the current family state, use that value for all
        currently_dc  = formula in st.session_state.double_charge_families
        changed_dc_rows = fam_edited[fam_edited["2+"] != currently_dc]
        if not changed_dc_rows.empty:
            new_dc = bool(changed_dc_rows.iloc[0]["2+"])
            if new_dc:
                st.session_state.double_charge_families.add(formula)
            else:
                st.session_state.double_charge_families.discard(formula)
            changed = True

        # Delete — tick ANY row → mark whole family; untick ANY row → unmark whole family
        pending = st.session_state.get("_pending_delete", set())
        currently_pending = formula in pending
        any_ticked   = bool(fam_edited["Delete"].any())
        any_unticked = not bool(fam_edited["Delete"].all())

        if not currently_pending and any_ticked:
            # At least one box just got ticked — select the whole family
            if "_pending_delete" not in st.session_state:
                st.session_state._pending_delete = set()
            st.session_state._pending_delete.add(formula)
            changed = True
        elif currently_pending and any_unticked:
            # At least one box just got unticked — deselect the whole family
            st.session_state._pending_delete.discard(formula)
            changed = True

    if changed:
        st.rerun()

    # ── Delete button ─────────────────────────────────────────────────────────────
    pending = st.session_state.get("_pending_delete", set())
    if pending:
        families = ", ".join(sorted(pending))
        col_del, col_cancel = st.columns([1, 4])
        with col_del:
            if st.button(f"🗑 Delete {families}", type="primary"):
                for formula in pending:
                    st.session_state.formula_list = [
                        f for f in st.session_state.formula_list if f != formula
                    ]
                    st.session_state.scaling_factors.pop(formula, None)
                    st.session_state.dc_scaling_factors.pop(formula, None)
                    st.session_state.double_charge_families.discard(formula)
                st.session_state._pending_delete = set()
                st.rerun()
        with col_cancel:
            if st.button("Cancel"):
                st.session_state._pending_delete = set()
                st.rerun()

    # Rebuild with confirmed stable state
    iso_df = _build_iso_df(formulas)

    # ── Build plot ────────────────────────────────────────────────────────────────
    sigma         = st.session_state.peak_width
    fig           = go.Figure()
    shown_legends : set[str] = set()

    for formula in formulas:
        colour   = colour_map[formula]
        fam_rows = iso_df[iso_df["Formula"] == formula]
        dc       = formula in st.session_state.double_charge_families

        for _, row in fam_rows.iterrows():
            # ── 1+ peak (always shown) ────────────────────────────────────────
            show_leg = formula not in shown_legends
            if show_leg:
                shown_legends.add(formula)

            fig.add_trace(_make_gaussian_trace(
                centre      = row["m/z"],
                amplitude   = row["_scaled"],
                sigma       = sigma,
                colour      = colour,
                name        = formula,
                showlegend  = show_leg,
                legendgroup = formula,
                hover_text  = f"Label: {row['Label']}<br>Abundance: {row['Abundance (%)']:.3f}%",
            ))

            # ── 2+ peak (shown as dashed outline when 2+ is ticked) ──────────
            if dc:
                dc_name = formula + " (2+)"
                show_dc_leg = dc_name not in shown_legends
                if show_dc_leg:
                    shown_legends.add(dc_name)

                # Dashed, lighter version of the same colour
                half_width = 4.5 * sigma
                x_dc = np.linspace(row["m/z (2+)"] - half_width, row["m/z (2+)"] + half_width, _GAUSSIAN_POINTS)
                y_dc = _gaussian_y(x_dc, row["_scaled_dc"], row["m/z (2+)"], sigma)
                fig.add_trace(go.Scatter(
                    x=x_dc, y=y_dc,
                    mode="lines",
                    fill="tozeroy",
                    fillcolor=_hex_to_rgba(colour, 0.10),
                    line=dict(color=colour, width=1.5, dash="dash"),
                    name=dc_name,
                    legendgroup=dc_name,
                    showlegend=show_dc_leg,
                    hovertemplate=(
                        f"<b>{dc_name}</b><br>"
                        f"m/z: %{{x:.4f}}<br>"
                        f"Intensity: %{{y:.4f}}<br>"
                        f"Label: {row['Label']} (2+)<br>"
                        f"Abundance: {row['Abundance (%)']:.3f}%"
                        "<extra></extra>"
                    ),
                ))

    # Spectrum overlay
    if st.session_state.spectrum_data is not None:
        mz_arr, int_arr = st.session_state.spectrum_data
        fig.add_trace(go.Scatter(
            x=mz_arr, y=int_arr,
            mode="lines",
            name="Mass Spectrum",
            line=dict(color="black", width=1.2),
            opacity=0.85,
        ))

    fig.update_layout(
        xaxis_title="m/z",
        yaxis_title="Intensity (scaled)",
        hovermode="closest",
        autosize=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(t=60, b=60, l=60, r=20),
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
    )
    fig.update_xaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)
    fig.update_yaxes(showgrid=True, gridcolor="rgba(128,128,128,0.2)", zeroline=False, showline=True, linewidth=1, linecolor="gray", mirror=True)

    st.plotly_chart(fig)

elif st.session_state.active_tool == "m/z Lookup":
    st.subheader("Reverse Lookup Results")

    candidates = []

    if run_reverse and target_mz > 0:

        for species_group in all_species.values():
            candidates.extend(species_group)

        candidates = sorted(set(candidates))

    filtered_candidates = []
    
    if element_filter.strip():

        requested_elements = {
            e.strip().capitalize()
            for e in element_filter.split(",")
            if e.strip()
        }

        for formula in candidates:
            try:
                composition = Formula(formula).composition()

                species_elements = set(composition.keys())

                if requested_elements.intersection(species_elements):
                    filtered_candidates.append(formula)

            except Exception:
                continue

    candidates = filtered_candidates
    matches = []

    for form in candidates:
        try:
            df = isotope_rows(form)
            for _, row in df.iterrows():
                if abs(row["mz"] - target_mz) <= tolerance:
                    matches.append({
                        "Element": form,
                        "Isotope": row["label"],
                        "m/z": row["mz"],
                        "Abundance (%)": row["abundance"] * 100,
                        "Error": abs(row["mz"] - target_mz),
                    })
        except Exception:
            continue

    if matches:
        match_df = pd.DataFrame(matches).sort_values("Error")
        st.success(f"Found {len(match_df)} matches")
        st.dataframe(match_df, width="stretch")
    else:
        st.warning("No matches found within tolerance.")

elif st.session_state.active_tool == "Oxide Converter":
    st.header("Oxide → Element Converter Results")

    if uploaded_file is not None:

        try:
            df = pd.read_csv(uploaded_file)

            st.subheader("Input")
            st.dataframe(df, width="stretch")

            # Remove empty rows
            df = df.dropna(subset=[df.columns[0], df.columns[1]])

            elemental_totals = defaultdict(float)

            for _, row in df.iterrows():

                formula_str = str(row.iloc[0]).strip()
                wt_pct = float(row.iloc[1])

                formula = Formula(formula_str)

                total_mass = formula.mass

                for element, item in formula.composition().items():

                    elemental_totals[element] += (
                        wt_pct * item.mass / total_mass
                    )

            output_df = pd.DataFrame({
                "Element": list(elemental_totals.keys()),
                "WeightPercent": list(elemental_totals.values())
            })

            output_df.sort_values("Element", inplace=True)

            st.subheader("Elemental Composition")
            st.dataframe(
                output_df.round(6),
                width="stretch"
            )

            csv_export = output_df.to_csv(index=False)

            st.download_button(
                label="Download Element Standard CSV",
                data=csv_export,
                file_name="elemental_standard.csv",
                mime="text/csv",
            )

        except Exception as e:
            st.error(f"Conversion failed: {e}")
