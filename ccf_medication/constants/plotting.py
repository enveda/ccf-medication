"""
Plotting constants
"""
REMISSION_COLOR = "#7A9B76"  # green
ACTIVE_COLOR = "#048BA8"  # blue
BASE_COLOR = "#60548b"  # purple

PLOTTING_DRUG_FAMILIES_MAP = {
    'TNF-alpha inhibitors': r'TNF-$\alpha$ inhibitors',
    'Integrin-receptor antagonists': 'Integrin Inhibitors',
    'IL-12/23 inhibitor': 'IL-12/23 Inhibitors',
    'JAK inhibitors': 'JAK Inhibitors',
    'Immunomodulating antimetabolites': 'Immunomodulators',
    'Aminosalicylates (5-ASA drugs)': '5-ASA Drugs',
    'No Medication': 'No Medication',
}

MUTED_MODERN_ALL_DRUGS = [
    "#60548b",  # muted deep purple
    "#f3b55a",  # soft golden yellow
    "#c2445b",  # dusty raspberry
    "#ed7c32",  # muted orange
    "#b4b6cf",  # periwinkle blue
    "#9b8691",  # taupe gray
    "#000000",  # black
]

MUTED_MODERN_REM_VS_ACT = [
    "#7A9B76",  # mantis green
    "#048BA8",  # almost turquoise blue
] 