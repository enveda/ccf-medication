# Adjusted p-value threshold
ADJ_PVAL_THRESH = 0.05

# Fold change threshold
TX_FC_THRESH = 1 # for transcriptomics
PX_FC_THRESH = 0.5 # for proteomics

# Minimum number of patients and samples for the Med vs No Med model (Remission only patients)
MIN_PATIENTS_PER_DRUG_FAMILY = 10

# Minimum number of patients and samples for the Rem vs Active model (Remission and Active patients) 
# this means there will be at least 10 total patients in the rem vs active analysis
MIN_PATIENTS_PER_SEVERITY_GROUP = 5

# Scientific notation threshold (for plotting)
SCIENTIFIC_NOTATION_THRESHOLD = 3
