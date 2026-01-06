#Load_environment
source ~/.bashrc
conda activate blast_env

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#Constants
file_path = 'merged_vfdb_blast_results.xlsx'
max_components = 10  #maximum_PCA_components

#Species_keyword
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

#Load_merged_Excel_file
xls = pd.ExcelFile(file_path)

#presence/absence_matrix
all_sseq = set()
strain_map = {}

for sheet in xls.sheet_names:

    df_sheet = pd.read_excel(xls, sheet_name=sheet)
    sseqs = df_sheet['sseqid'].dropna().unique()

    strain_map[sheet] = set(sseqs)
    all_sseq.update(sseqs)
all_sseq = sorted(all_sseq)
pa = pd.DataFrame(0, index=xls.sheet_names, columns=all_sseq)

for sheet, sseqs in strain_map.items():
    pa.loc[sheet, list(sseqs)] = 1

#Extract_species_label

def get_label(sheet_name):

    if sheet_name.startswith('D__'):

        for sp in dickeya_species:

            if f'D__{sp}' in sheet_name:

                return 'D_' + sp

    if sheet_name.startswith('P__'):

        for sp in pecto_species:

            if f'P__{sp}' in sheet_name:

                return 'P_' + sp

    return None



pa['Label'] = pa.index.map(get_label)
pa = pa[pa['Label'].notna()].copy()

#List_of_labels_for_Dickeya_and_Pectobacterium_separately

d_labels = ['D_' + sp for sp in dickeya_species]
p_labels = ['P_' + sp for sp in pecto_species]

#Filter_to_present

d_labels = [lbl for lbl in d_labels if lbl in pa['Label'].values]
p_labels = [lbl for lbl in p_labels if lbl in pa['Label'].values]

labels = d_labels + p_labels

#Color_map
total_labels = labels
palette = sns.color_palette('hls', len(total_labels))
color_map = {lbl: palette[i] for i, lbl in enumerate(total_labels)}
palette = sns.color_palette('hls', len(labels))
color_map = {lbl: palette[i] for i, lbl in enumerate(labels)}


#PCA_function
def pca_for_label(label):

    strains = pa[pa['Label'] == label].index.tolist()

    n_samples = len(strains)

    if n_samples < 2:

        print(f"Skipping {label}: only {n_samples} strain(s)")

        return


    #Data_submatrix

    X = pa.loc[strains, all_sseq]
    Xs = StandardScaler().fit_transform(X)


    #Determine_PCA_components

    n_comp = min(max_components, n_samples, Xs.shape[1])

    if n_comp < 2:

        print(f"Skipping {label}: insufficient samples/features for PCA")

        return

    pca = PCA(n_components=n_comp)
    Xp = pca.fit_transform(Xs)
    var_pct = pca.explained_variance_ratio_ * 100
    pcs = [f'PC{i+1}' for i in range(n_comp)]



    #Feature_contributions
    R = np.abs(pca.components_.T)
    feat_cont = R / R.sum(axis=0)
    feat_df = pd.DataFrame(feat_cont, index=all_sseq, columns=pcs)
    ev = R.dot(pca.explained_variance_ratio_)
    feat_df['Explained_by_feature'] = ev / ev.sum() * var_pct.sum()
    feat_df.to_excel(f'featcont_{label}.xlsx', float_format='%.5f')

    print(f"Wrote featcont_{label}.xlsx")



    #2D_PCA_scatter_with_strain_labels
    coords = pd.DataFrame(Xp, index=strains, columns=pcs)
    plt.figure(figsize=(8,6))
    for st in strains:

        x, y = coords.loc[st, ['PC1', 'PC2']]

        plt.scatter(x, y, color=color_map[label], s=50, alpha=0.7)
        plt.text(x + 0.01, y + 0.01, st, fontsize=4, color='black')
    plt.xlabel(f'PC1 ({var_pct[0]:.1f}%)')
    plt.ylabel(f'PC2 ({var_pct[1]:.1f}%)')
    plt.title(f'2D PCA {label} (strain labels)')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'2D_{label}.png', dpi=1200)
    plt.close()
  
    print(f"Saved 2D_{label}.png")



        #Pairwise_PC1-4_plot 
    n_pair = min(4, n_comp)

    if len(strains) >= 3 and n_pair >= 2:

        rename_map = {pcs[i]: f"{pcs[i]} ({var_pct[i]:.1f}%)" for i in range(n_pair)}
        df_pair = coords.rename(columns=rename_map)

        #KDE_on_diagonal_to_avoid_histogram_bin_issues

        sns.pairplot(

            df_pair,
            vars=list(rename_map.values()),
            diag_kind='kde',
            corner=True,
            plot_kws={'alpha':0.7, 's':40}
        )

        plt.suptitle(f'Pairwise PCs (1–{n_pair}) {label}', y=1.02)
        plt.tight_layout()
        plt.savefig(f'pairplot_{label}.png', dpi=1200)
        plt.close()

        print(f"Saved pairplot_{label}.png")

    else:

        print(f"Skipping pairplot for {label}: too few strains ({len(strains)}) or PCs ({n_pair})")



#PCA_for_each_label
for lbl in labels:

    pca_for_label(lbl)



print('All species PCAs complete.')('All species PCAs complete.')

