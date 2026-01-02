#!/bin/bash
#SBATCH --job-name=WG_speciesPCA
#SBATCH --output=WG_speciesPCA.out
#SBATCH --error=WG_speciesPCA.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G

# Author: Kalyani Bhandari
# Project: Soft Rot Pectobacteriaceae comparative genomics
# Purpose: Whole-genome protein PCA at species level across Dickeya, Pectobacterium, and Erwinia
# HPC: Tempest (SLURM)

echo "Starting Whole Genome Species-level PCA job"

# Load conda and activate environment
source ~/.bashrc
conda activate blast_env

# Move to working directory
cd /home/g89x126/Labwork_Bioinformatics || exit 1

# Create output folder
mkdir -p WG_PCA_output
cd WG_PCA_output || exit 1

#PYTHON CODE STARTS HERE#
python << 'EOF'
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

# =======================
# 0) Paths
# =======================
DATA_DIR = "/home/g89x126/Labwork_Bioinformatics/Protein_faa"

# =======================
# 1) FASTA file detection
# =======================
faa_files = [f for f in os.listdir(DATA_DIR) if f.endswith((".faa", ".fa"))]
print("Total protein files found:", len(faa_files))

# =======================
# 2) Extract gene names
# =======================
def extract_gene_name(header):
    header = header.replace(">", "").strip()
    before_bracket = header.split("[")[0].strip()
    tokens = before_bracket.split()
    return tokens[-1].strip(",.;:()[]") if tokens else None

# =======================
# 3) Presence/absence matrix
# =======================
strain_to_genes = {}
all_genes = set()

for fname in faa_files:
    strain = fname.replace(".faa", "").replace(".fa", "")
    genes = set()
    with open(os.path.join(DATA_DIR, fname)) as f:
        for line in f:
            if line.startswith(">"):
                g = extract_gene_name(line)
                if g:
                    genes.add(g)
                    all_genes.add(g)
    strain_to_genes[strain] = genes

all_genes = sorted(all_genes)
PA = pd.DataFrame(0, index=strain_to_genes.keys(), columns=all_genes)
for strain, genes in strain_to_genes.items():
    PA.loc[strain, list(genes)] = 1

print("Total unique gene names:", len(all_genes))

# =======================
# 4) Species definitions
# =======================
dickeya_species = [
    'aquatica','chrysanthemi','dadantii','dianthicola','fangzhongdai',
    'lacustris','oryzae','parazeae','poaceiphila','solani','zea'
]

pecto_species = [
    'actinidae','aquaticum','araliae','aroidearum','betavasculorum',
    'brasiliense','carotovorum','colocasium','fontis','jejuense',
    'odoriferum','parmentieri','parvum','peruviense','polaris',
    'polonicum','punjabense','quasiaquaticum','versatile',
    'wasabiae','zantedeschiae'
]

def classify_strain(name):
    if name.startswith("D__"):
        for sp in dickeya_species:
            if f"D__{sp}" in name:
                return ("Dickeya", "D_" + sp)
    if name.startswith("P__"):
        for sp in pecto_species:
            if f"P__{sp}" in name:
                return ("Pectobacterium", "P_" + sp)
    if name.startswith("Erwinia_"):
        return ("Erwinia", "E_Erwinia")
    return (None, None)

PA["Genus"], PA["SpeciesLabel"] = zip(*PA.index.map(classify_strain))
PA = PA[PA["Genus"].notna()].copy()

# =======================
# 5) PCA subsets
# =======================
subsets = {
    "Dickeya": PA[PA["Genus"] == "Dickeya"].index.tolist(),
    "Pectobacterium": PA[PA["Genus"] == "Pectobacterium"].index.tolist(),
    "Erwinia": PA[PA["Genus"] == "Erwinia"].index.tolist(),
    "All": PA.index.tolist(),
}

unique_species = sorted(PA["SpeciesLabel"].unique())
palette = sns.color_palette("hls", len(unique_species))
color_map = dict(zip(unique_species, palette))

# =======================
# 6) PCA function
# =======================
def run_pca(name, strains):
    X = PA.loc[strains, all_genes].values
    if X.shape[0] < 2:
        return

    X_std = StandardScaler().fit_transform(X)
    pca = PCA(n_components=min(10, *X.shape))
    Xp = pca.fit_transform(X_std)
    var = pca.explained_variance_ratio_ * 100

    coords = pd.DataFrame(Xp[:, :2], index=strains, columns=["PC1","PC2"])
    coords["SpeciesLabel"] = PA.loc[strains, "SpeciesLabel"]
    coords.to_csv(f"WG_PCAcoords_{name}.csv")

    plt.figure(figsize=(8,6))
    for sp in unique_species:
        sub = coords[coords["SpeciesLabel"] == sp]
        if not sub.empty:
            plt.scatter(sub.PC1, sub.PC2, s=40, alpha=0.7,
                        color=color_map[sp], label=sp)
    plt.xlabel(f"PC1 ({var[0]:.1f}%)")
    plt.ylabel(f"PC2 ({var[1]:.1f}%)")
    plt.title(f"PCA – {name}")
    plt.legend(fontsize=6, bbox_to_anchor=(1.04,1))
    plt.grid(True)
    plt.savefig(f"WG_PCA_{name}.png", dpi=1200)
    plt.close()

for g in subsets:
    run_pca(g, subsets[g])

print("ALL PCA RUN COMPLETE")
EOF
#PYTHON_CODE_ENDS_HERE

echo "Whole Genome Species-level PCA job finished."
