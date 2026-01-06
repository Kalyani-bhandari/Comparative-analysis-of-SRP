#!/bin/bash
#SBATCH --job-name=WG_PFAM_HMM_PCA
#SBATCH --output=WG_PFAM_HMM_PCA.out
#SBATCH --error=WG_PFAM_HMM_PCA.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G

echo "Starting Whole Genome Pfam HMM"

#Load_environment
source ~/.bashrc
conda activate blast_env

#work_directory
cd /home/username/Labwork_Bioinformatics

PFAM_DIR="Pfam_HMM"
HMM_RES="HMM_results"
PROT_DIR="Protein_faa"

mkdir -p "$PFAM_DIR" "$HMM_RES"

echo "Using protein directory: $PROT_DIR"
echo "Using Pfam directory: $PFAM_DIR"
echo "Results in: $HMM_RES"

#1_Download_Pfam-A_and_hmmpress

cd "$PFAM_DIR"

if [[ ! -f "Pfam-A.hmm" ]]; then
    echo "Downloading Pfam-A HMM…"
    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
else
    echo "Pfam-A.hmm already exists. Skipping download."
fi

cd /home/username/Labwork_Bioinformatics

#2_Run_hmmscan_for_every_.faa_file

echo "Running hmmscan on all protein files"

FILES=( ${PROT_DIR}/*.faa ${PROT_DIR}/*.fa )

echo "Found ${#FILES[@]} protein files."

for FILE in "${FILES[@]}"; do

    STRAIN=$(basename "$FILE")
    STRAIN="${STRAIN%.*}"

    OUT="${HMM_RES}/${STRAIN}.domtblout"

    if [[ -s "$OUT" ]]; then
        echo "Skipping $STRAIN (already processed)"
        continue
    fi

    echo "Running hmmscan for $STRAIN …"

    hmmscan --cpu ${SLURM_CPUS_PER_TASK} \
            --domtblout "$OUT" \
            "${PFAM_DIR}/Pfam-A.hmm" \
            "$FILE" \
            > "${HMM_RES}/${STRAIN}.log"
done

echo "hmmscan complete for ALL strains"
