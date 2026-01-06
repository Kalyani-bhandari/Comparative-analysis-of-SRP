#!/bin/bash
#
#SBATCH --account=priority-janakjoshi
#SBATCH --partition=priority
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=14-00:00:00
#SBATCH --job-name=Orthofinder_vfdb
#SBATCH --output=Orthofinder_vfdb-%j.out
#SBATCH --error=Orthofinder_vfdb-%j.err
#SBATCH --mail-user=
#SBATCH --mail-type=ALL

#Author: Kalyani Bhandari
#Project: Soft Rot Pectobacteriaceae comparative genomics
#Purpose: BLASTP-based screening of all strain proteomes against VFDB
#Database: Virulence Factor Database (VFDB)
#HPC: Tempest (SLURM)

#Load_conda
source ~/.bashrc
conda activate blast_env

#Load_BLAST+
module load BLAST+/2.14.0-gompi-2022b

#Directories
input_dir=~/Labwork_Bioinformatics/Protein_faa
vfdb_fasta="$input_dir/vfdb_protein.fasta"
blast_db="$input_dir/vfdb_protein_db"
output_dir=~/Labwork_Bioinformatics/BLAST_results_vfdb
mkdir -p "$output_dir"

#Make_BLAST_database
makeblastdb -in "$vfdb_fasta" -dbtype prot -out "$blast_db"

#Filters
e_value_cutoff="1e-3"
perc_identity_cutoff="40"

#Run_BLAST
for fasta_file in "$input_dir"/*.faa; do
    [[ -f "$fasta_file" ]] || continue

    fasta_name=$(basename "$fasta_file" .faa)
    result_dir="$output_dir/$fasta_name"
    mkdir -p "$result_dir"

    echo "Running BLASTP on: $fasta_name"

    blastp -query "$fasta_file" -db "$blast_db" \
        -out "$result_dir/${fasta_name}_results.txt" \
        -outfmt 6 \
        -evalue "$e_value_cutoff" \
        -num_threads 32

    #Header_for_filtered_results
    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" \
        > "$result_dir/${fasta_name}_filtered_results.txt"

    #Identity_and_e-value_filters
    awk -v pid="$perc_identity_cutoff" -v ev="$e_value_cutoff" \
        '$3 >= pid && $11 <= ev' \
        "$result_dir/${fasta_name}_results.txt" \
        >> "$result_dir/${fasta_name}_filtered_results.txt"

    echo "Finished processing: $fasta_name"
done

echo "All VFDB BLAST analyses completed."
