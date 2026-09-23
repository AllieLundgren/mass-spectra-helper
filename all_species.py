# Mass spectrometry species database for laser ablation / laser ionization (LA-LI)
# For metals, critical minerals, rare earths, battery materials
#
# Hand-curated lists come first; the rule-generated families at the bottom add
# every plausible cluster, oxide, hydride, halide, etc. that LA-LI plasmas form.
# Every entry is a plain formula string that molmass can parse. Duplicates
# (same elemental composition, e.g. "Fe(OH)2" and "FeO2H2") are fine — the app
# keeps the first name it sees, so the curated names win.

from itertools import combinations
from math import ceil

elements = [
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
    "Au","Hg","Tl","Pb","Bi",
    "Th","U",
]

rare_earths = [
    "Sc","Y","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu",
]

# Small clusters and dimers
dimers_trimers = [
    # Carbon
    "C2", "C3", "C4", "C5",
    # Silicon
    "Si2", "Si3",
    # Oxygen
    "O2", "O3",
    # Nitrogen
    "N2",
    # Metal dimers (common in sputtering/evaporation)
    "Fe2", "Cu2", "Ni2", "Al2", "Mg2", "Ti2", "Co2", "Zn2",
    "Ag2", "Au2", "Pt2", "Pd2", "Mo2", "W2", "Ta2",
    # Rare earth dimers
    "Ce2", "Nd2", "Dy2", "Yb2", "La2",
    # Mixed trimers
    "Fe3", "Cu3", "Al3",
]

# Monoxides
monoxides = [
    # Transition metals
    "SiO","TiO","VO","CrO","MnO","FeO","CoO","NiO","CuO","ZnO",
    "MoO","WO","ReO","OsO","IrO","PtO","RuO","RhO","PdO",
    # Main group
    "AlO","MgO","CaO","SrO","BaO","ScO","YO",
    "ZrO","NbO","HfO","TaO","GeO","SnO","PbO",
    # Lanthanides
    "LaO","CeO","PrO","NdO","SmO","EuO","GdO","TbO","DyO",
    "HoO","ErO","TmO","YbO","LuO",
    # Actinides
    "ThO","UO",
    # Alkali/alkaline earth
    "K2O","Na2O","Li2O","Rb2O","Cs2O",
]

# Metal dioxides
dioxides = [
    "SiO2","TiO2","ZrO2","HfO2","CeO2","AlO2","MoO2","WO2",
    "GeO2","SnO2","PbO2","MnO2","FeO2","CoO2","NiO2","CuO2","ZnO2",
    "NbO2","TaO2","ReO2","UO2","ThO2","VO2",
]

# Sesquioxides and other mixed oxides
mixed_oxides = [
    "Fe2O3","Al2O3","Cr2O3","Ga2O3","La2O3","Ce2O3","Nd2O3",
    "Sm2O3","Eu2O3","Gd2O3","Dy2O3","Ho2O3","Yb2O3","Y2O3",
    "Sc2O3","Ti2O3","V2O3","Mn2O3","Fe2O4","Cu2O3",
    "P2O5","P4O10","As2O5","Sb2O5","B2O3",
    "Mn3O4","Fe3O4","Co3O4",
    "MoO3","WO3","UO3","U3O8","V2O5","Nb2O5","Ta2O5",
]

# Nitrides
nitrides = [
    "BN","AlN","GaN","InN","ZnN",
    "TiN","VN","CrN","MnN","FeN","CoN","NiN","CuN",
    "ZrN","NbN","MoN","TaN","WN","HfN",
    "Si3N4","Al2N3","Mg3N2",
]

# Carbides
carbides = [
    "SiC","TiC","VC","CrC","MoC","WC","ZrC","NbC","TaC","HfC",
    "Fe3C","Fe4C","Co2C","Ni2C",
    "B4C","B2C",  # Boron carbide fragments
    "LiC","Li2C",
]

# Sulfides, selenides, tellurides
chalcogenides = [
    "ZnS","CdS","HgS","FeS","FeS2","CoS","NiS","CuS","MoS2","WS2",
    "As2S3","Sb2S3",
    "ZnSe","CdSe","HgSe",
    "ZnTe","CdTe","HgTe",
    "PbS","PbSe","PbTe",
    "SnS","SnS2","SnSe",
]

# Phosphides, arsenides, antimonides
pnictides = [
    "AlP","GaP","InP","ZnP",
    "AlAs","GaAs","InAs",
    "AlSb","GaSb","InSb",
    "Mg2P","Ca2P","Zn2P",
]

# Metal silicates and aluminosilicates (partial formulas that appear in MS)
silicates_aluminates = [
    "SiOH","SiO3","SiO4","Si2O5","Si3O6","Si4O8",
    "AlOH","AlO3","AlO4",
    "MgSiO3","CaSiO3","BaSiO3","ZnSiO3",
    "Al2SiO5",
]

# Lithium compounds (battery relevance)
lithium_compounds = [
    "LiH","LiO","LiO2","Li2O","Li2O2",
    "LiC","Li2C","Li2C2",
    "LiN","Li3N",
    "LiS","Li2S","Li2S2",
    "LiF","LiCl","LiBr","LiI",
    "LiOH",
    # Layered oxides and phosphates
    "LiCoO2","LiFeO2","LiNiO2","LiMnO2",
    "LiFePO4","LiMnPO4","LiCoPO4",
    "LiSiO2","LiAlO2","LiZrO2",
    "LiTiO2","LiVO2","LiCrO2",
    "LiMn2O4","Li2CO3","LiPF6","LiPO2F2",
]

# Sodium compounds (battery materials)
sodium_compounds = [
    "NaO","Na2O","Na2O2",
    "NaC","Na2C",
    "NaS","Na2S",
    "NaF","NaCl","NaBr",
    "NaOH",
    "NaCoO2","NaFeO2","NaNiO2","NaMnO2",
    "NaFePO4","NaMnPO4",
]

# Metal alloys and intermetallics (especially relevant for battery anodes/cathodes)
alloys_intermetallics = [
    # Lithium alloys
    "LiSi","Li2Si","Li4Si","Li5Si",
    "LiAl","Li2Al","Li3Al",
    "LiMg","Li2Mg",
    "LiNi","Li2Ni","Li3Ni",
    "LiCo","Li2Co","Li3Co",
    "LiZn","Li2Zn","Li4Zn",
    "LiFe","Li2Fe",
    "LiSn","Li2Sn","Li4Sn",
    "LiPb","Li2Pb",
    # Sodium alloys
    "NaSi","Na2Si",
    "NaAl","Na2Al",
    "NaSn","Na2Sn",
    # Magnesium alloys
    "MgSi","Mg2Si",
    "MgZn","Mg2Zn",
    "MgCu","Mg2Cu",
    "MgNi","Mg2Ni",
    # Copper alloys
    "CuZn","Cu2Zn","Cu3Zn",
    "CuNi","Cu3Ni",
    "CuAl","Cu3Al","Cu9Al4",
    # Other common intermetallics
    "FeAl","Fe2Al","Fe3Al",
    "NiAl","Ni3Al",
    "TiAl","Ti3Al",
    "CoAl","Co2Al",
    "MnAl","Mn3Al",
    "FeNi","Fe3Ni","Fe4Ni",
]

# Ternary and quaternary oxides (relevant to battery materials)
complex_oxides = [
    # Spinel-like
    "MgAl2O4","ZnAl2O4","CoAl2O4",
    "MgFe2O4","ZnFe2O4","CoFe2O4",
    # Perovskite-like fragments
    "LaAlO3","LaTiO3","LaFeO3","LaCoO3","LaNiO3",
    "CaTiO3","SrTiO3","BaTiO3",
    # Pyrochlore-like
    "La2Zr2O7","La2Ti2O7",
    # Garnet-like
    "Y3Al5O12",  # may appear as fragments
    "Li7La3Zr2O12",  # LLZO solid electrolyte
    # NASICON-like
    "LiTi2P3O12","LiZr2P3O12",
    # Olivine fragments (already covered in LiFePO4, etc.)
]

# Phosphates
phosphates = [
    "PO","PO2","PO3","PO4",
    "P2O5","P3O5","P4O10",
    "HPO3","H2PO4","HPO4",
    "AlPO4","FePO4","CoPO4","NiPO4","ZnPO4",
    "LaPO4","CePO4",
]

# Hydroxides
hydroxides = [
    "OH","H2O",
    "LiOH","NaOH","KOH",
    "Al(OH)","Al(OH)3",
    "Fe(OH)2","Fe(OH)3","FeOOH",
    "Co(OH)2","Ni(OH)2","Cu(OH)2","Zn(OH)2",
    "Mg(OH)2","Ca(OH)2","Ba(OH)2",
    "Si(OH)4",
    "Ti(OH)4",
]

# Hydrides
hydrides = [
    "LiH","NaH","KH","CaH2","MgH2",
    "AlH3",
    "SiH4",
    "PH3",
    "AsH3",
    "B2H6",
]

# Fluorides and halides
halides = [
    "LiF","LiCl","LiBr","LiI",
    "NaF","NaCl","NaBr","NaI",
    "KF","KCl","KBr","KI",
    "MgF2","MgCl2",
    "AlF3","AlCl3",
    "ZnF2","ZnCl2",
    "CaF2","CaCl2",
    "FeF2","FeF3","FeCl2","FeCl3",
    "CoF2","CoCl2",
    "NiF2","NiCl2",
    "CuF2","CuCl2",
    "TiF4","TiCl4",
    "ZrF4","ZrCl4",
]

# Rare earth specific compounds (fluorides, chlorides, sesquioxides, phosphates)
rare_earth_compounds = (
    [f"{re_}F3" for re_ in rare_earths]
    + [f"{re_}Cl3" for re_ in rare_earths]
    + [f"{re_}2O3" for re_ in rare_earths]
    + [f"{re_}PO4" for re_ in rare_earths]
)

# Nitrogen oxides
nitrogen_oxides = [
    "NO","NO2","NO3","N2O","N2O3","N2O4","N2O5",
    "N2","N3",
]

# Fragment ions and common MS artifacts (residual gas, surface contamination)
common_fragments = [
    "H","H2","H3",
    "OH","H2O","H3O",
    "O","O2","O3",
    "N","N2","N3",
    "NH","NH2","NH3","NH4",
    "C","C2","C3","CO","CO2","HCO","HCO2",
    "CN","NCO","HCN",
    "Si","SiO","SiO2",
    "P","PO","PO2",
    "S","SO","SO2","SO3","SO4","HS","H2S",
    "Se","SeO","SeO2",
    "Cl","ClO","ClO2","HCl",
    "F","HF","Br","BrO","HBr","I","IO","HI",
]


# ═════════════════════════════════════════════════════════════════════════════
# Rule-generated families
# ═════════════════════════════════════════════════════════════════════════════

NONMETALS = {"H","He","C","N","O","F","Ne","P","S","Cl","Ar","Se","Br","Kr","I","Xe"}
HALOGENS  = ["F","Cl","Br","I"]
ALKALI    = ["Li","Na","K","Rb","Cs"]

# Metals + metalloids (B, Si, Ge, As, Sb, Te) — everything that forms M_x X_y ions
METALS = [e for e in elements if e not in NONMETALS]

# Highest common oxidation state — bounds how many O / halogen atoms a metal
# can reasonably carry (a little over this is allowed for sub-/super-oxide ions)
MAX_OX = {
    "Li":1, "Be":2, "B":3,  "Na":1, "Mg":2, "Al":3, "Si":4, "K":1,  "Ca":2,
    "Sc":3, "Ti":4, "V":5,  "Cr":6, "Mn":7, "Fe":3, "Co":3, "Ni":3, "Cu":2,
    "Zn":2, "Ga":3, "Ge":4, "As":5, "Rb":1, "Sr":2, "Y":3,  "Zr":4, "Nb":5,
    "Mo":6, "Ru":8, "Rh":4, "Pd":4, "Ag":2, "Cd":2, "In":3, "Sn":4, "Sb":5,
    "Te":6, "Cs":1, "Ba":2, "La":3, "Ce":4, "Pr":4, "Nd":3, "Sm":3, "Eu":3,
    "Gd":3, "Tb":4, "Dy":3, "Ho":3, "Er":3, "Tm":3, "Yb":3, "Lu":3, "Hf":4,
    "Ta":5, "W":6,  "Re":7, "Os":8, "Ir":6, "Pt":6, "Au":3, "Hg":2, "Tl":3,
    "Pb":4, "Bi":5, "Th":4, "U":6,
}


def _f(*parts):
    """Build a formula string from (symbol, count) pairs, e.g. ("Na", 2), ("Cl", 1) → "Na2Cl"."""
    return "".join(f"{sym}{n if n > 1 else ''}" for sym, n in parts if n > 0)


def _ox(m):
    return MAX_OX.get(m, 3)


# Homonuclear metal clusters M2–M5
metal_clusters = [_f((m, n)) for m in METALS for n in range(2, 6)]

# Nonmetal / metalloid clusters beyond what metal_clusters covers
nonmetal_clusters = (
    [_f((e, n)) for e in ("Si", "B") for n in range(6, 9)]
    + [_f((e, n)) for e in ("P", "S") for n in range(2, 9)]
    + [_f(("Se", n)) for n in range(2, 7)]
    + [_f((e, n)) for e in ("H", "N", "O") for n in range(2, 5)]
)

# Carbon clusters, hydrocarbons and small C–N / C–O fragments (organics, binders, carbon coatings)
carbon_clusters_hydrocarbons = (
    [_f(("C", n)) for n in range(2, 21)]
    + [_f(("C", n), ("H", h)) for n in range(1, 13) for h in range(1, 2 * n + 3)]
    + [_f(("C", n), ("N", k)) for n in range(1, 7) for k in (1, 2)]
    + [_f(("C", n), ("O", k)) for n in range(1, 5) for k in (1, 2)]
)

# Oxides M_xO_y, x = 1–3, y up to ~half the max oxidation state per metal + 1
metal_oxides = [
    _f((m, x), ("O", y))
    for m in METALS for x in (1, 2, 3)
    for y in range(1, ceil(x * _ox(m) / 2) + 2)
]

# Hydrides MH_n, M2H, oxy-hydrides / hydroxides MO_yH_n, M2OH, alkali (MOH)nM series
metal_hydrides_hydroxides = (
    [_f((m, 1), ("H", k)) for m in METALS for k in range(1, min(3, _ox(m) + 1) + 1)]
    + [_f((m, 2), ("H", 1)) for m in METALS]
    + [
        _f((m, 1), ("O", y), ("H", k))
        for m in METALS for y in range(1, ceil(_ox(m) / 2) + 2) for k in (1, 2)
    ]
    + [_f((m, 2), ("O", 1), ("H", 1)) for m in METALS]
    + [_f((a, n + 1), ("O", n), ("H", n)) for a in ALKALI for n in range(1, 5)]
)

# Nitrides M_xN_y, x, y = 1–2
metal_nitrides = [_f((m, x), ("N", y)) for m in METALS for x in (1, 2) for y in (1, 2)]

# Carbides MC1–4, M2C1–3
metal_carbides = (
    [_f((m, 1), ("C", y)) for m in METALS for y in range(1, 5)]
    + [_f((m, 2), ("C", y)) for m in METALS for y in range(1, 4)]
)

# Halides M_xX_y, x = 1–3, plus alkali-halide cluster series (MX)n and (MX)nM
# — e.g. Na2Cl = (NaCl)Na, the classic NaCl cluster ion
metal_halides = [
    _f((m, x), (X, y))
    for X in HALOGENS for m in METALS for x in (1, 2, 3)
    for y in range(1, min(x * _ox(m) + 1, 8) + 1)
] + [
    _f((a, n + extra), (X, n))
    for a in ALKALI for X in HALOGENS for n in range(1, 9) for extra in (0, 1)
]

# Oxyhalides MOX, MOX2, MO2X (e.g. LaOCl, ZrOCl2)
oxyhalides = [
    _f((m, 1), ("O", y), (X, k))
    for X in HALOGENS for m in METALS for y, k in ((1, 1), (1, 2), (2, 1))
]

# Sulfides / selenides / tellurides M_xZ_y, x = 1–2, y = 1–3
metal_chalcogenides = [
    _f((m, x), (Z, y))
    for Z in ("S", "Se", "Te") for m in METALS if m != Z
    for x in (1, 2) for y in (1, 2, 3)
]

# Phosphides / arsenides / antimonides M_xPn_y, x, y = 1–2
metal_pnictides = [
    _f((m, x), (Pn, y))
    for Pn in ("P", "As", "Sb") for m in METALS if m != Pn
    for x in (1, 2) for y in (1, 2)
]

# Metal + oxyanion fragments: borates, carbonates, nitrates, silicates, phosphates, sulfates
# M_kZO_y, k = 1–2, y = 1–4 (e.g. LiPO3, Na2SO4, CaCO3, Na2NO3 = (NaNO3)Na)
oxyanion_salts = [
    _f((m, k), (Z, 1), ("O", y))
    for Z in ("B", "C", "N", "Si", "P", "S") for m in METALS if m != Z
    for k in (1, 2) for y in range(1, 5)
]

# Heteronuclear metal clusters AB, A2B, AB2 for every metal pair (alloys, intermetallics)
mixed_metal_clusters = [
    _f((a, i), (b, j))
    for a, b in combinations(METALS, 2) for i, j in ((1, 1), (2, 1), (1, 2))
]

# Mixed-metal oxides ABO1–4 for every metal pair (e.g. LiCoO2, LaAlO3, MgSiO3)
mixed_metal_oxides = [
    _f((a, 1), (b, 1), ("O", y))
    for a, b in combinations(METALS, 2) for y in range(1, 5)
]

# Argon adducts — only relevant if argon is used as a buffer/carrier gas.
# Off by default in the app's lookup (see DEFAULT_OFF_GROUPS).
argon_adducts = [_f((m, 1), ("Ar", 1)) for m in METALS] + [
    "Ar2","ArH","ArH2","Ar2H","ArO","ArOH","ArN","ArC","ArN2","ArO2",
]


# Compile all into master dictionary (curated first so their names win on duplicates)
all_species = {
    "elements": elements,
    "dimers_trimers": dimers_trimers,
    "monoxides": monoxides,
    "dioxides": dioxides,
    "mixed_oxides": mixed_oxides,
    "nitrides": nitrides,
    "carbides": carbides,
    "chalcogenides": chalcogenides,
    "pnictides": pnictides,
    "silicates_aluminates": silicates_aluminates,
    "lithium_compounds": lithium_compounds,
    "sodium_compounds": sodium_compounds,
    "alloys_intermetallics": alloys_intermetallics,
    "complex_oxides": complex_oxides,
    "phosphates": phosphates,
    "hydroxides": hydroxides,
    "hydrides": hydrides,
    "halides": halides,
    "rare_earth_compounds": rare_earth_compounds,
    "nitrogen_oxides": nitrogen_oxides,
    "common_fragments": common_fragments,
    # ── rule-generated ──
    "metal_clusters": metal_clusters,
    "nonmetal_clusters": nonmetal_clusters,
    "carbon_clusters_hydrocarbons": carbon_clusters_hydrocarbons,
    "metal_oxides": metal_oxides,
    "metal_hydrides_hydroxides": metal_hydrides_hydroxides,
    "metal_nitrides": metal_nitrides,
    "metal_carbides": metal_carbides,
    "metal_halides": metal_halides,
    "oxyhalides": oxyhalides,
    "metal_chalcogenides": metal_chalcogenides,
    "metal_pnictides": metal_pnictides,
    "oxyanion_salts": oxyanion_salts,
    "mixed_metal_clusters": mixed_metal_clusters,
    "mixed_metal_oxides": mixed_metal_oxides,
    "argon_adducts": argon_adducts,
}

# Groups the m/z lookup leaves unticked until the user turns them on
DEFAULT_OFF_GROUPS = {"argon_adducts"}
