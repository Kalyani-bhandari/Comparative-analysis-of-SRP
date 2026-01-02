#!/bin/bash
#SBATCH --job-name=iuc_operon_blast
#SBATCH --output=iuc_operon_blast_%j.out
#SBATCH --error=iuc_operon_blast_%j.err
#SBATCH --ntasks=1
#SBATCH --partition=nextgen
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

# Author: Kalyani Bhandari
# Project: Soft Rot Pectobacteriaceae comparative genomics
# Purpose: BLASTP-based identification of aerobactin (iuc) operon genes
#          across Dickeya and Pectobacterium type strains
# HPC: Tempest (SLURM)

# Activate conda environment
source ~/.bashrc
conda activate blast_env

# -------------------- Paths --------------------
BASE_DIR="/home/g89x126/Labwork_Bioinformatics/iuc_operon"
QUERY="$BASE_DIR/iuc_operon.faa"

DICT_DIR="$BASE_DIR/Dickeya_typestrain"
PECTO_DIR="$BASE_DIR/Pectobacterium_typestrain"

OUT_DIR="$BASE_DIR/blast_results"
mkdir -p "$OUT_DIR"

# Output files
DIC_OUT="$OUT_DIR/Dickeya_iuc_operon_all.tsv"
PEC_OUT="$OUT_DIR/Pectobacterium_iuc_operon_all.tsv"

# Reset output files
> "$DIC_OUT"
> "$PEC_OUT"

# -------------------- Parameters --------------------
EVALUE=1e-5

# Header
HEADER="strain\tqseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqcovs\tannotation"
echo -e "$HEADER" > "$DIC_OUT"
echo -e "$HEADER" > "$PEC_OUT"

# -------------------- Dickeya --------------------
echo "Running BLAST for Dickeya type strains..."
for SUBJECT in "$DICT_DIR"/*.faa; do
    BASENAME=$(basename "$SUBJECT" .faa)
    echo "Processing $BASENAME ..."
    
    blastp -query "$QUERY" -subject "$SUBJECT" \
        -evalue $EVALUE \
        -num_threads 8 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
        | awk -v strain="$BASENAME" 'BEGIN{FS=OFS="\t"} {print strain,$0}' >> "$DIC_OUT"
done

# -------------------- Pectobacterium --------------------
echo "Running BLAST for Pectobacterium type strains..."
for SUBJECT in "$PECTO_DIR"/*.faa; do
    BASENAME=$(basename "$SUBJECT" .faa)
    echo "Processing $BASENAME ..."
    
    blastp -query "$QUERY" -subject "$SUBJECT" \
        -evalue $EVALUE \
        -num_threads 8 \
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
        | awk -v strain="$BASENAME" 'BEGIN{FS=OFS="\t"} {print strain,$0}' >> "$PEC_OUT"
done

echo "All BLAST analyses completed."
echo "Results saved to:"
echo "  $DIC_OUT"
echo "  $PEC_OUT"
