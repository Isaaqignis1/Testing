"""shared genus colours and plot ordering, used by every figure script."""

# archaeal genera — colour by clade
color_map = {
    # Haloarchaea — blues
    "Haloarcula":           "#08306b",
    "Halobacterium":        "#2171b5",
    "Haloferax":            "#6baed6",
    "Halorubrum":           "#c0dbf5",
    # Methanogens — greens
    "Methanobacterium_D":   "#00b92e",
    "Methanobrevibacter_A": "#73ff00",
    "Methanococcus":        "#001f0d",
    "Methanoculleus":       "#006161",
    "Methanohalophilus":    "#6de781",
    "Methanothermobacter":  "#c7ff98",
    "Methanosphaera":       "#62976a",
    # Ammonia oxidisers — purples
    "Nitrosopelagicus":     "#4a148c",
    "Nitrososphaera":       "#ce93d8",
    # Thermoacidophiles — warm
    "Metallosphaera":       "#bd9a00",
    "Saccharolobus":        "#ff9d5c",
    "Sulfolobus":           "#d32f2f",
    "Sulfuracidifex":       "#a623e2",
    # Thermococci — reds
    "Thermococcus":         "#580101",
    "Thermococcus_B":       "#B48D99",
}

genus_order = [
    "Haloarcula", "Halobacterium", "Haloferax", "Halorubrum",
    "Methanobacterium_D", "Methanobrevibacter_A", "Methanoculleus",
    "Methanothermobacter", "Methanosphaera", "Methanococcus",
    "Methanohalophilus",
    "Nitrosopelagicus", "Nitrososphaera",
    "Metallosphaera", "Saccharolobus", "Sulfolobus", "Sulfuracidifex",
    "Thermococcus", "Thermococcus_B",
]

# bacterial control colours
control_color_map = {
    "Campylobacter": "#0033dd",
    "Escherichia":   "#888888",
    "Brucella":      "#5e2e8a",
}

control_order = ["Brucella", "Escherichia", "Campylobacter"]
