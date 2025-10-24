# ------------------------------------------------------------
#### CONSTANTS ####
# ------------------------------------------------------------

BASE_METADATA_COLUMNS_OF_INTEREST = [
    "patient_id",
    "diagnosis",
    "simple_tissue",
    "medication",
    "endo_category",
    "macroscopic_appearance",
    "birth_year",
    "age_at_sample_collection",
    "moa",
    "sex",
]

CATEGORIES_OF_INTEREST = {"endo_category": ["remission", "mild", "moderate", "severe"]}

METADATA_COLUMNS_OF_INTEREST = list(
    set(BASE_METADATA_COLUMNS_OF_INTEREST + list(CATEGORIES_OF_INTEREST.keys()))
)

METADATA_2024_COLS_RENAME_MAP = {
    "MACROSCOPIC_APPEARANCE": "macroscopic_appearance",
    "DEIDENTIFIED_MASTER_PATIENT_ID": "patient_id",
    "SAMPLE_COLLECTED_DATE": "sample_collected_date",
    "SEX": "sex",
    "MEDICATION_AT_INDEX": "raw_medication",
    "DIAGNOSIS": "diagnosis",
    "tissue": "simple_tissue",
    "ENDO_CATEGORY": "endo_category",
    "sample_id": "transcriptomics",
    "BIRTH_YEAR": "birth_year",
    "MOA": "moa",
}


DRUG_FAMILIES = {
    "Integrin-receptor antagonists": ["Vedolizumab", "Natalizumab"],
    "TNF-alpha inhibitors": [
        "Infliximab (Remicade)",
        "Infliximab (Inflectra)",
        "Infliximab (Renflexis)",
        "Infliximab (Avsola)",
        "Adalimumab (Humira)",
        "Certolizumab Pegol",
        "Golimumab",
    ],
    "IL-12/23 inhibitor": ["Ustekinumab"],
    "JAK inhibitors": ["Tofacitinib", "Upadacitinib"],
    "Immunomodulating antimetabolites": [
        "Azathioprine",
        "Mercaptopurine",
        "Methotrexate",
    ],
    "Calcineurin inhibitors": ["Tacrolimus", "Cyclosporine"],
    "Aminosalicylates (5-ASA drugs)": ["Mesalamine", "Balsalazide", "Sulfasalazine"],
    "Combination Therapy": ["Combination Therapy"],
    "No Medication": ["No Medication"],
}

NO_MEDICATION = "No Medication"
PATIENT_ID_COL = "patient_id"
SAMPLE_ID_COL = "sample_id"
DRUG_FAMILY_COL = "drug_family"
DEFAULT_EXCLUDE_FAMILIES = ["Combination Therapy"]

P_VALUE_COLUMN = "adjusted_p_value"
FOLD_CHANGE_COLUMN = "coefficient"
FEATURE_COLUMN = "feature"

LOGGER_NAME = "ccf_medication"

PVAL_COL_REM_VS_ACT = "adj_p_val_rem_vs_active"
FC_COL_REM_VS_ACT = "coef_rem_vs_active"


PVAL_COL_MED_VS_NO_MED = "adj_p_val_med_vs_no_med"
FC_COL_MED_VS_NO_MED = "coef_med_vs_no_med"
