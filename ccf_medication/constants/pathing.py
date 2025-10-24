import os
from pathlib import Path

# ------------------------------------------------------------
#### DATA PATHS ####
# ------------------------------------------------------------
HOME_FOLDER = str(Path(__file__).resolve().parent.parent.parent)
## Original DataFrame paths

## TODO: These paths need to be added
METADATA_PRE_2024_PATH = None
METADATA_2024_PATH = None
TX_PATH = None
PX_PATH = None

## Root data directory
DATA_DIR = os.path.join(HOME_FOLDER, "data")

## Processed DataFrame paths
PROCESSED_DATA_DIR = os.path.join(DATA_DIR, "processed_dataframes")

TX_DATAFRAME_PATH = os.path.join(
    PROCESSED_DATA_DIR, "transcriptomics_dataframe.parquet"
)

PX_DATAFRAME_PATH = os.path.join(PROCESSED_DATA_DIR, "proteomics_dataframe.parquet")

REMISSION_TX_DATAFRAME_PATH = os.path.join(
    PROCESSED_DATA_DIR, "transcriptomics_dataframe_remission.parquet"
)

ACTIVE_TX_DATAFRAME_PATH = os.path.join(
    PROCESSED_DATA_DIR, "transcriptomics_dataframe_active.parquet"
)

REMISSION_PX_DATAFRAME_PATH = os.path.join(
    PROCESSED_DATA_DIR, "proteomics_dataframe_remission.parquet"
)

ACTIVE_PX_DATAFRAME_PATH = os.path.join(
    PROCESSED_DATA_DIR, "proteomics_dataframe_active.parquet"
)

## Gene column paths
TX_GENE_COLS = os.path.join(PROCESSED_DATA_DIR, "transcriptomics_gene_cols.txt")

PX_PROTEIN_COLS = os.path.join(PROCESSED_DATA_DIR, "proteomics_protein_cols.txt")

## Patient table paths
REMISSION_PATIENT_TABLE_PATH = os.path.join(
    PROCESSED_DATA_DIR, "remission_patient_table.parquet"
)
ACTIVE_PATIENT_TABLE_PATH = os.path.join(
    PROCESSED_DATA_DIR, "active_patient_table.parquet"
)

## Results paths
RESULTS_DIR = os.path.join(DATA_DIR, "results")

MED_VS_NO_MED_RESULTS_DIR = os.path.join(RESULTS_DIR, "med_vs_no_med")

MED_VS_NO_MED_TX_DIR = os.path.join(MED_VS_NO_MED_RESULTS_DIR, "transcriptomics")

MED_VS_NO_MED_PX_DIR = os.path.join(MED_VS_NO_MED_RESULTS_DIR, "proteomics")

AGG_MED_VS_NO_MED_TX_RESULTS_PATH = os.path.join(
    MED_VS_NO_MED_TX_DIR, "agg_med_vs_no_med_transcriptomics_results.parquet"
)

AGG_MED_VS_NO_MED_PX_RESULTS_PATH = os.path.join(
    MED_VS_NO_MED_PX_DIR, "agg_med_vs_no_med_proteomics_results.parquet"
)

## Disease severity results paths
REM_VS_ACT_RESULTS_DIR = os.path.join(RESULTS_DIR, "rem_vs_active")

REM_VS_ACT_TX_DIR = os.path.join(REM_VS_ACT_RESULTS_DIR, "transcriptomics")

REM_VS_ACT_PX_DIR = os.path.join(REM_VS_ACT_RESULTS_DIR, "proteomics")

AGG_REM_VS_ACT_TX_RESULTS_PATH = os.path.join(
    REM_VS_ACT_TX_DIR, "agg_rem_vs_active_transcriptomics_results.parquet"
)

AGG_REM_VS_ACT_PX_RESULTS_PATH = os.path.join(
    REM_VS_ACT_PX_DIR, "agg_rem_vs_active_proteomics_results.parquet"
)

## Drug vs no med table Path
MED_VS_NO_MED_SIGNIF_TABLE_PATH = os.path.join(
    PROCESSED_DATA_DIR, "med_vs_no_med_signif_table.parquet"
)

## Disease severity table Path
REM_VS_ACT_SIGNIF_TABLE_PATH = os.path.join(
    PROCESSED_DATA_DIR, "rem_vs_active_signif_table.parquet"
)

## Responder gene strip plots
RESPONEDER_ANALYSIS_DIR = os.path.join(RESULTS_DIR, "responder_analysis")
TX_RESPONSE_ANALYSIS_DIR = os.path.join(RESPONEDER_ANALYSIS_DIR, "transcriptomics")
PX_RESPONSE_ANALYSIS_DIR = os.path.join(RESPONEDER_ANALYSIS_DIR, "proteomics")
RESPONSE_TABLE_PATH = os.path.join(RESPONEDER_ANALYSIS_DIR, "responder_table.parquet")

# Create directories
for path in [
    DATA_DIR,
    PROCESSED_DATA_DIR,
    RESULTS_DIR,
    MED_VS_NO_MED_RESULTS_DIR,
    MED_VS_NO_MED_TX_DIR,
    MED_VS_NO_MED_PX_DIR,
    RESPONEDER_ANALYSIS_DIR,
    TX_RESPONSE_ANALYSIS_DIR,
    PX_RESPONSE_ANALYSIS_DIR,
]:
    os.makedirs(path, exist_ok=True)
