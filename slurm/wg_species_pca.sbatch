#!/bin/bash
#SBATCH --job-name=WG_speciesPCA
#SBATCH --output=WG_speciesPCA.out
#SBATCH --error=WG_speciesPCA.err
#SBATCH --time=3-00:00:00
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G

#Author: Kalyani Bhandari
#Project: Soft Rot Pectobacteriaceae comparative genomics
#Purpose: Whole-genome protein PCA at species level across Dickeya, Pectobacterium, and Erwinia
#HPC: Tempest (SLURM)

echo "Starting Whole Genome Species-level PCA job"

#activate_environment
source ~/.bashrc
conda activate blast_env

#working_directory
cd /home/....../.......

#output_folder
mkdir -p WG_PCA_output
cd WG_PCA_output

#PYTHON_CODE_STARTS_HERE#

python << 'EOF'
import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#0_Paths
DATA_DIR = "/home/g89x126/Labwork_Bioinformatics/Protein_faa"

#1_FASTA_file_detection
faa_files = [
    f for f in os.listdir(DATA_DIR)
    if f.endswith(".faa") or f.endswith(".fa")
]

print("Total protein files found:", len(faa_files))

#2_Extract_gene_name_from_FASTA_header
def extract_gene_name(header):
    header = header.replace(">", "").strip()
    before_bracket = header.split("[")[0].strip()
    tokens = before_bracket.split()

    if len(tokens) == 0:
        return None

    #Prefer_predicted_gene_names_first
    last = tokens[-1].strip(",.;:()[]")
    return last if last else None

#3_Build_presence/absence_matrix
strain_to_genes = {}
all_genes = set()

for fname in faa_files:
    strain = fname.replace(".faa", "").replace(".fa", "")
    full_path = os.path.join(DATA_DIR, fname)

    genes = set()
    with open(full_path) as f:
        for line in f:
            if line.startswith(">"):
                gene_name = extract_gene_name(line)
                if gene_name:
                    genes.add(gene_name)
                    all_genes.add(gene_name)

    strain_to_genes[strain] = genes

all_genes = sorted(all_genes)
print("Total unique gene names:", len(all_genes))

PA = pd.DataFrame(0, index=strain_to_genes.keys(), columns=all_genes)
for strain, genes in strain_to_genes.items():
    PA.loc[strain, list(genes)] = 1

#4_Species_definitions
dickeya_species = [
    'aquatica','chrysanthemi','dadantii','dianthicola',
    'fangzhongdai','lacustris','oryzae','parazeae',
    'poaceiphila','solani','zea'
]

pecto_species = [
    'actinidae','aquaticum','araliae','aroidearum',
    'betavasculorum','brasiliense','carotovorum',
    'colocasium','fontis','jejuense','odoriferum',
    'parmentieri','parvum','peruviense','polaris',
    'polonicum','punjabense','quasiaquaticum',
    'versatile','wasabiae','zantedeschiae'
]

#5_Species_extraction_INCLUDING_ERWINIA
def classify_strain(name):
    #Dickeya
    if name.startswith("D__"):
        for sp in dickeya_species:
            if f"D__{sp}" in name:
                return ("Dickeya", "D_" + sp)

    #Pectobacterium
    if name.startswith("P__"):
        for sp in pecto_species:
            if f"P__{sp}" in name:
                return ("Pectobacterium", "P_" + sp)

    #Erwinia_generic
    if name.startswith("Erwinia_"):
        return ("Erwinia", "E_Erwinia")

    #fallback
    return (None, None)

PA["Genus"], PA["SpeciesLabel"] = zip(*PA.index.map(classify_strain))
PA = PA[PA["Genus"].notna()].copy()

print("Strains with assigned genus:", PA.shape[0])

#6_Subsets
subsets = {
    "Dickeya": PA[PA["Genus"] == "Dickeya"].index.tolist(),
    "Pectobacterium": PA[PA["Genus"] == "Pectobacterium"].index.tolist(),
    "Erwinia": PA[PA["Genus"] == "Erwinia"].index.tolist(),
    "All": PA.index.tolist(),
}

#7_Color_map
unique_species = sorted(PA["SpeciesLabel"].unique())
palette = sns.color_palette("hls", len(unique_species))
color_map = {sp: palette[i] for i, sp in enumerate(unique_species)}

#8_PCA_Function
def run_pca(name, strains):
    print(f"\nRunning PCA: {name} (n={len(strains)})")

    X = PA.loc[strains, all_genes].values
    n_samples, n_features = X.shape

    if n_samples < 2 or n_features < 2:
        print("Skipping (not enough data)")
        return

    #Standardize
    X_std = StandardScaler().fit_transform(X)

    n_comp = min(10, n_samples, n_features)
    pca = PCA(n_components=n_comp)
    Xp = pca.fit_transform(X_std)
    var = pca.explained_variance_ratio_ * 100
    pcs = [f"PC{i+1}" for i in range(n_comp)]

    coords = pd.DataFrame(Xp, index=strains, columns=pcs)
    coords["SpeciesLabel"] = PA.loc[strains, "SpeciesLabel"]

    #Save_coordinates
    coords.to_csv(f"WG_PCAcoords_{name}.csv")

    #Scree_plot
    plt.figure(figsize=(6,4))
    plt.plot(range(1,n_comp+1), var, "o-")
    plt.xlabel("PC")
    plt.ylabel("Variance (%)")
    plt.title(f"Scree Plot – {name}")
    plt.grid(True)
    plt.savefig(f"WG_Scree_{name}.png", dpi=900)
    plt.close()

    #PCA_Plot
    plt.figure(figsize=(8,6))
    for sp in unique_species:
        sub = coords[coords["SpeciesLabel"] == sp]
        if not sub.empty:
            plt.scatter(sub["PC1"], sub["PC2"], s=40, alpha=.7,
                        color=color_map[sp], label=sp)

    plt.xlabel(f"PC1 ({var[0]:.1f}%)")
    plt.ylabel(f"PC2 ({var[1]:.1f}%)")
    plt.title(f"PCA – {name}")
    plt.legend(fontsize=6, bbox_to_anchor=(1.04,1))
    plt.grid(True)
    plt.savefig(f"WG_PCA_{name}.png", dpi=1200)
    plt.close()

    #Feature_contributions
    R = np.abs(pca.components_.T)
    feature_cont = R / R.sum(axis=0)

    feat_df = pd.DataFrame(feature_cont, index=all_genes, columns=pcs)
    ev = R.dot(pca.explained_variance_ratio_)
    feat_df["ExplainedVariance_Total"] = ev / ev.sum() * var.sum()
    feat_df.to_excel(f"WG_FeatureContrib_{name}.xlsx")

    print(f"Finished {name} | PC1 variance = {var[0]:.1f}%")

#9_RUN_ALL_SUBSETS

for group in ["Dickeya", "Pectobacterium", "Erwinia", "All"]:
    run_pca(group, subsets[group])

print("\nALL PCA RUN COMPLETE ✓")
EOF
#PYTHON_CODE_ENDS_HERE#

echo "Whole Genome Species-level PCA job finished."

