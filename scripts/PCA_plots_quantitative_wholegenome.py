#Load_environment
source ~/.bashrc
conda activate blast_env

#run_the_script(python name.py)
import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from mpl_toolkits.mplot3d import Axes3D

#DIRECTORIES
dom_dir = "HMM_results"
pfam_annotation_file = "Pfam_HMM/Pfam-A.clans.tsv" 
out_dir = "PCA_quantitative_results"
os.makedirs(out_dir, exist_ok=True)

#SPECIES_GROUPING
dickeya_species = [
    'aquatica','chrysanthemi','dadantii','dianthicola',
    'fangzhongdai','lacustris','oryzae','parazeae',
    'poaceiphila','solani','zea','undicola'
]

pecto_species = [
    'actinidae','aquaticum','araliae','aroidearum',
    'betavasculorum','brasiliense','carotovorum',
    'colocasium','fontis','jejuense','odoriferum',
    'parmentieri','parvum','peruviense','polaris',
    'polonicum','punjabense','quasiaquaticum',
    'versatile','wasabiae','zantedeschiae','atrosepticum'
]

#1_LOAD_PFAM_ANNOTATION_FILE
pfam_annot = pd.read_csv(
    pfam_annotation_file,
    sep="\t",
    header=None,
    names=["Acc", "Clan", "ID", "Desc"]
)

pfam_annot["Acc"] = pfam_annot["Acc"].str.strip()
pfam_annot = pfam_annot.drop_duplicates("Acc").set_index("Acc")

print("Pfam annotation entries loaded:", len(pfam_annot))

#2_PARSE_HMMER_domtblout

def parse_domtblout_quantitative(filepath):
    """
    Extracts domain -> highest bitscore for each strain.
    """
    domain_scores = {}

    with open(filepath, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue

            cols = line.strip().split()
            if len(cols) < 13:
                continue

            domain = cols[1].split(".")[0]
            bitscore = float(cols[11])  # correct name

            # Keep the maximum bit score
            if domain not in domain_scores:
                domain_scores[domain] = bitscore
            else:
                domain_scores[domain] = max(domain_scores[domain], bitscore)

    return domain_scores

#3_LOAD_ALL_domtblout_FILES

dom_files = sorted([f for f in os.listdir(dom_dir) if f.endswith(".domtblout")])
print("Found", len(dom_files), "HMM result files.")

strain_scores = {}
all_domains = set()

for domfile in dom_files:
    strain = domfile.replace(".domtblout", "")
    path = os.path.join(dom_dir, domfile)

    score_dict = parse_domtblout_quantitative(path)
    strain_scores[strain] = score_dict
    all_domains.update(score_dict.keys())

all_domains = sorted(all_domains)
print("Total unique domains with quantitative scores:", len(all_domains))

#4_BUILD_SCORE_MATRIX

df = pd.DataFrame(0.0, index=strain_scores.keys(), columns=all_domains)

for strain, doms in strain_scores.items():
    for d, score in doms.items():
        df.loc[strain, d] = score

#5_SPECIES_LABELING

def detect_species(s):
    if s.startswith("D__"):
        for sp in dickeya_species:
            if f"D__{sp}" in s:
                return "D_" + sp

    if s.startswith("P__"):
        for sp in pecto_species:
            if f"P__{sp}" in s:
                return "P_" + sp

    if s.startswith("Erwinia_"):
        m = re.match(r"Erwinia_([A-Za-z]+)", s)
        return "E_" + m.group(1) if m else "E_sp"

    return "Unknown"

df["SpeciesLabel"] = df.index.map(detect_species)
df = df[df["SpeciesLabel"] != "Unknown"]  # remove unknown strains

feature_cols = [c for c in df.columns if c != "SpeciesLabel"]

#6_PCA_FUNCTION

def run_pca(name, subset):
    print("\n=== Running Quantitative PCA:", name, "===")

    data = df.loc[subset, feature_cols]
    if data.shape[0] < 3:
        print("Too few samples for PCA.")
        return

    X = StandardScaler().fit_transform(data)
    n_comp = min(10, data.shape[0], data.shape[1])

    pca = PCA(n_components=n_comp)
    PC = pca.fit_transform(X)

    pcs = [f"PC{i+1}" for i in range(n_comp)]

    coords = pd.DataFrame(PC, index=subset, columns=pcs)
    coords["SpeciesLabel"] = df.loc[subset, "SpeciesLabel"]
    coords.to_csv(f"{out_dir}/PCA_quantitative_coords_{name}.csv")

    #Scree_Plot
    var = pca.explained_variance_ratio_ * 100
    plt.figure()
    plt.plot(range(1, n_comp+1), var, "o-")
    plt.ylabel("% Variance Explained")
    plt.xlabel("PC")
    plt.title(f"Scree Plot - {name}")
    plt.savefig(f"{out_dir}/Scree_{name}.png", dpi=900, bbox_inches="tight")
    plt.close()

    #Feature_Contributions
    R = np.abs(pca.components_.T)
    feat_df = pd.DataFrame(R, index=feature_cols, columns=pcs)
    feat_df["TotalContribution"] = R.dot(pca.explained_variance_ratio_)

    #Merge_Pfam_annotation
    merged = feat_df.merge(
        pfam_annot,
        left_index=True,
        right_index=True,
        how="left"
    )

    merged.to_excel(f"{out_dir}/FeatureContrib_{name}_ANNOTATED.xlsx")
    merged["TotalContribution"].sort_values(ascending=False).to_csv(
        f"{out_dir}/TopDomains_{name}_ANNOTATED.csv"
    )

    #Color_Map
    species_labels = coords["SpeciesLabel"].unique()
    colors = plt.cm.tab20(np.linspace(0, 1, len(species_labels)))

    color_map = {}
    for sp, col in zip(species_labels, colors):
        if sp.startswith("E_"):
            color_map[sp] = (0.2, 0.2, 0.2)  # grey for all Erwinias
        else:
            color_map[sp] = col

    #2D_PCA_PLOT
    plt.figure(figsize=(9,7))

    for sp in species_labels:
        sub = coords[coords["SpeciesLabel"] == sp]

        # OPTIONAL drop shadow
        plt.scatter(
            sub["PC1"] + 0.15,
            sub["PC2"] - 0.15,
            s=100,
            color="grey",
            alpha=0.15
        )

        # REAL glossy bubble
        plt.scatter(
            sub["PC1"], sub["PC2"],
            s=100,
            color=color_map[sp],
            alpha=0.70,
            linewidth=0.7,
            marker="o",
            label=sp
        )

    plt.title(f"PCA Quantitative - {name}")
    plt.xlabel(f"PC1 ({var[0]:.1f}%)")
    plt.ylabel(f"PC2 ({var[1]:.1f}%)")
    plt.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.02,0.5))
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.savefig(f"{out_dir}/PCA_quantitative_{name}.png", dpi=900, bbox_inches="tight")
    plt.close()

    #LABELED_PCA
    plt.figure(figsize=(10,8))

    for sp in species_labels:
        sub = coords[coords["SpeciesLabel"] == sp]
        plt.scatter(
            sub["PC1"], sub["PC2"],
            s=100,
            color=color_map[sp],
            alpha=0.70,
            edgecolor="black",
            linewidth=0.7
        )

        for idx, row in sub.iterrows():
            plt.text(row["PC1"], row["PC2"], idx, fontsize=5)

    plt.title(f"PCA Quantitative (Labeled) - {name}")
    plt.grid(True, linestyle="--", alpha=0.3)
    plt.savefig(f"{out_dir}/PCA_quantitative_{name}_labeled.png", dpi=900, bbox_inches="tight")
    plt.close()


    #3D_PCA_PLOT
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111, projection="3d")

    for sp in species_labels:
        sub = coords[coords["SpeciesLabel"] == sp]

        ax.scatter(
            sub["PC1"], sub["PC2"], sub["PC3"],
            s=120,
            color=color_map[sp],
            alpha=0.75,
            edgecolor="black",
            linewidth=0.7
        )

    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_zlabel(f"PC3 ({var[2]:.1f}%)")
    ax.set_title(f"3D PCA Quantitative - {name}")

    plt.savefig(f"{out_dir}/PCA_quantitative_{name}_3D.png", dpi=900, bbox_inches="tight")
    plt.close()
  
#RUN_PCA
allsp = df.index.tolist()
dickeya = df[df["SpeciesLabel"].str.startswith("D_")].index.tolist()
pects = df[df["SpeciesLabel"].str.startswith("P_")].index.tolist()

run_pca("All", allsp)
run_pca("Dickeya", dickeya)
run_pca("Pectobacterium", pects)

print("\n Quantitative PCA complete. Output in PCA_quantitative_results/")
