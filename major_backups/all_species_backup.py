# Expanded mass spectrometry species database
# For metals, critical minerals, rare earths, battery materials

elements = [
    "H","He","Li","Be","B","C","N","O","F","Ne",
    "Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
    "Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Ru","Rh","Pd","Ag","Cd",
    "In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy",
    "Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt",
    "Au","Hg","Tl","Pb","Bi"
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

# Monoxides (your original list expanded)
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
    # Alkali/alkaline earth
    "K2O","Na2O","Li2O","Rb2O","Cs2O",
]

# Dioxides (your original list expanded)
dioxides = [
    "SiO2","TiO2","ZrO2","HfO2","CeO2","AlO2","MoO2","WO2","CO2",
    "GeO2","SnO2","PbO2","MnO2","FeO2","CoO2","NiO2","CuO2","ZnO2",
    "NbO2","TaO2","ReO2","UO2","ThO2","VO2",
]

# Sesquioxides and other mixed oxides
mixed_oxides = [
    "Fe2O3","Al2O3","Cr2O3","Ga2O3","La2O3","Ce2O3","Nd2O3",
    "Sm2O3","Eu2O3","Gd2O3","Dy2O3","Ho2O3","Yb2O3","Y2O3",
    "Sc2O3","Ti2O3","V2O3","Mn2O3","Fe2O4","Cu2O3",
    "P2O5","P4O10","As2O5","Sb2O5","B2O3",
    "Mn3O4","Fe3O4",
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
    "B4C","B2C","B","C3",  # Boron carbide fragments
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
    "SiOH","SiO3","SiO4",
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

# Rare earth specific compounds
rare_earth_compounds = [
    # Oxides, dioxides already covered above
    # Fluorides
    "LaF3","CeF3","PrF3","NdF3","SmF3","EuF3","GdF3","DyF3","YbF3","LuF3",
    # Chlorides
    "LaCl3","CeCl3","PrCl3","NdCl3","SmCl3","EuCl3","GdCl3","DyCl3","YbCl3","LuCl3",
    # Mixed oxides
    "RE2O3",  # General rare earth oxide
    # Phosphates
    "LaPO4","CePO4","PrPO4","NdPO4","SmPO4","EuPO4","GdPO4","DyPO4","YbPO4","LuPO4",
]

# Nitrogen oxides
nitrogen_oxides = [
    "NO","NO2","N2O","N2O3","N2O4","N2O5",
    "N2","N3",
]

# Fragment ions and common MS artifacts
common_fragments = [
    "H","H2","H3",
    "OH","H2O","H3O",
    "O","O2","O3",
    "N","N2","N3",
    "C","C2","C3","CO","CO2",
    "Si","SiO","SiO2",
    "P","PO","PO2",
    "S","SO","SO2","SO3",
    "Cl","ClO",
    "CN","NCO",
]

# Compile all into master dictionary
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
}

