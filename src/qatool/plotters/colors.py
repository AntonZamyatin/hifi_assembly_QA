"""Shared color palettes and default tolerance constants for the report plots."""

PID_TIE_EPS_DEFAULT = 0.0001
BEST_PID_DELTA_DEFAULT = 0.001

SUPPORT_COLORS = {
    "asm1_only": "#b2182b",   # strong red
    "asm1_best": "#f4a3c2",   # light pink
    "asm2_only": "#1a9641",   # strong green
    "asm2_best": "#a6d96a",   # light green
    "common":    "#a3a395",   # grey-green
    "other":     "#424241",   # dark grey
}

# Occurrence classification colors: indices 0-4 = 1/2/3/4/5+ best loci (non-common),
# indices 5-9 = same counts for common reads (pale versions).
OCCURRENCE_COLORS = [
    "#2166ac",  # 1 locus
    "#92c5de",  # 1 locus, common
    "#1a9641",  # 2 loci
    "#a6d96a",  # 2 loci, common
    "#d7191c",  # 3 loci
    "#f4a582",  # 3 loci, common
    "#e6ab02",  # 4 loci
    "#fee08b",  # 4 loci, common
    "#252525",  # 5+ loci
    "#969696",  # 5+ loci, common
]
