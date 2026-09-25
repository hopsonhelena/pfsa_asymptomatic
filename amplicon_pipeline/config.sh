# config.sh: settings shared by every pipeline step.
# Edit the lines marked EDIT for your system. Scripts source this file from the
# pipeline directory (`source ./config.sh`), so run commands from there.

# ---- Paths (EDIT) ----
export PROJECT_DIR="$PWD"                                          
export DATA_DIR="${PROJECT_DIR}/data/runs"   # folder containing the sequencing run folders; can point elsewhere

# Joint reference: Ensembl GRCh38 primary assembly +
# P. falciparum 3D7, PlasmoDB release 62.
# See README.md#reference-genome for sources and preparation instructions.
# Contig names must match data/bed/
export REF="EDIT: path to concatenated reference FASTA"

export BED_DIR="${PROJECT_DIR}/data/bed"
export OUTPUT_DIR="${PROJECT_DIR}/output"
export WORK_DIR="${OUTPUT_DIR}/runs"
export HUMAN_DIR="${OUTPUT_DIR}/human"
export PARASITE_DIR="${OUTPUT_DIR}/parasite"
export COVERAGE_TABLE="${OUTPUT_DIR}/all_runs_mean_coverage.txt"

# ---- Sample lists ----
export METADATA="${PROJECT_DIR}/data/sample_metadata.csv"
# Optional filter lists: delete the file or leave it empty to skip 
export NO_METADATA_LIST="${PROJECT_DIR}/data/lackmetadata_samples.csv"
export DUPLICATES_LIST="${PROJECT_DIR}/data/duplicate_samples.csv"
export FEBRILE_LIST="${PROJECT_DIR}/data/febrile_samples.csv"

# ---- Sequencing runs (EDIT) ----
# One folder per run under ${DATA_DIR}. FASTQs are named ${run%R}X*_R{1,2}_001.fastq.gz
# (for example run 24118R contains 24118X100_..._R1_001.fastq.gz).
ALL_RUNS=(24118R 24404R 24448R 24466R 24494R 24538R 25579R 25631R)

# SLURM options (EDIT)
SLURM_ACCOUNT="EDIT: your Slurm account"
SLURM_PARTITION="EDIT: your Slurm partition"
N_RUNS=${#ALL_RUNS[@]}  # number of sequencing runs (sets the slurm array range: 0 to N_RUNS-1)


# Input FASTQs read from ${DATA_DIR}/${run}/${FASTQ_SUBDIR}
FASTQ_SUBDIR="Fastq"

# Derived files are written to ${WORK_DIR}/${run}/...
TRIMMED_SUBDIR="trimmed.fastq"
FASTQC_SUBDIR="fastqc"
BAM_SUBDIR="bams"
DEPTH_SUBDIR="mpileup_depth"

# ---- Amplicons (names of the bed files in data/bed, without .bed) ----
PF_AMPS=(acs8_1 acs8_2 ama1 sera2 trap phos csp)
HUMAN_AMPS=(hbs ab odel atp2b4 g6pd_202 g6pd_376)

# ---- Trimming (01_trim.sh, cutadapt) ----
ADAPTER_R1="CTGTCTCTTATACACATCT"
ADAPTER_R2="CTGTCTCTTATACACATCT"
MIN_LENGTH=20
