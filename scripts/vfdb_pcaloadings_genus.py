
#Load_environment
source ~/.bashrc
conda activate blast_env

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#Load_merged_Excel

file_path = 'merged_vfdb_blast_results.xlsx'
xls = pd.ExcelFile(file_path)

#Build_presence/absence_matrix_based_on_sseqid

all_sseqids = set()

strain_sseq_dict = {}

for sheet in xls.sheet_names:

    df_sheet = pd.read_excel(xls, sheet_name=sheet)
    sseqs = df_sheet['sseqid'].dropna().unique()
    strain_sseq_dict[sheet] = set(sseqs)
    all_sseqids.update(sseqs)
all_sseqids = sorted(all_sseqids)
presence_absence = pd.DataFrame(0, index=xls.sheet_names, columns=all_sseqids)

for sheet, sseqs in strain_sseq_dict.items():
    presence_absence.loc[sheet, list(sseqs)] = 1


#Define_species

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



#Extract_species_from_sheet_name

def extract_species(sheet):

    if sheet.startswith('D__'):
        for sp in dickeya_species:

            if f'D__{sp}' in sheet:
                return sp

    if sheet.startswith('P__'):
        for sp in pecto_species:

            if f'P__{sp}' in sheet:
                return sp

    return None



presence_absence['Species'] = presence_absence.index.map(extract_species)
df = presence_absence[presence_absence['Species'].notna()].copy()

#Prepare_subsets

subsets = {

    'Dickeya': df[df.index.str.startswith('D__')].index.tolist(),

    'Pectobacterium': df[df.index.str.startswith('P__')].index.tolist(),

    'Both': df.index.tolist()

}

#Assign_colors

species_list = sorted(df['Species'].unique())

species_label_list = [

    ('D_' if sp in dickeya_species else 'P_') + sp

    for sp in species_list

]

palette = sns.color_palette('hls', len(species_label_list))

label_color_map = {label: palette[i] for i, label in enumerate(species_label_list)}


#plotting_function

def species_pca(name, indices):

    #Subset_data

    sub_pa = df.loc[indices, all_sseqids]
    raw_species = df.loc[indices, 'Species']

    species_label = [

        ('D_' if sp in dickeya_species else 'P_') + sp

        for sp in raw_species

    ]


    #Standardize

    scaler = StandardScaler()
    X = scaler.fit_transform(sub_pa)

    #PCA_with_10_components

    n_comp = 10
    pca = PCA(n_components=n_comp)
    Xp = pca.fit_transform(X)
    var_pct = pca.explained_variance_ratio_ * 100

    #Build_coords_DataFrame

    pcs = [f'PC{i+1}' for i in range(n_comp)]
    coords_df = pd.DataFrame(Xp, index=indices, columns=pcs)
    coords_df['SpeciesLabel'] = species_label
    coords_df.to_csv(f'pca10_{name}.csv', index_label='Strain')

    #Scree_plot

    plt.figure(figsize=(6,4))
    plt.plot(range(1,n_comp+1), var_pct, 'o-', linewidth=1)
    plt.xlabel('PC')
    plt.ylabel('Variance Explained (%)')
    plt.title(f'Scree Plot ({name})')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'scree_{name}.png', dpi=1200)
    plt.close()


    #Cumulative_variance_plot

    cumvar = np.cumsum(var_pct)
    plt.figure(figsize=(6,4))
    plt.plot(range(1,n_comp+1), cumvar, 'o-', linewidth=1)
    for t in [50, 70]: plt.axhline(t, linestyle='--', label=f'{t}%')
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Variance Explained (%)')
    plt.title(f'Cumulative Variance ({name})')
    plt.legend(fontsize=6)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'cumvar_{name}.png', dpi=1200)
    plt.close()


    #Feature_contributions

    r = np.abs(pca.components_.T)
    feat_cont = r / r.sum(axis=0)
    feat_df = pd.DataFrame(feat_cont, index=all_sseqids, columns=pcs)
    ev = r.dot(pca.explained_variance_ratio_)
    feat_df['Explained_by_feature'] = ev/ev.sum() * var_pct.sum()
    feat_df.to_excel(f'featcont_{name}.xlsx', float_format='%.5f')


    #2D_PCA_no_labels

    plt.figure(figsize=(8,6))
    for label in sorted(species_label_list):

        sub = coords_df[coords_df['SpeciesLabel']==label]

        if not sub.empty:

            plt.scatter(sub['PC1'], sub['PC2'], color=label_color_map[label], label=label, s=50, alpha=0.7)

    plt.xlabel(f'PC1 ({var_pct[0]:.1f}%)')
    plt.ylabel(f'PC2 ({var_pct[1]:.1f}%)')
    plt.title(f'2D PCA {name} - No Labels')
    handles, labels = plt.gca().get_legend_handles_labels()
    pairs = sorted(zip(labels, handles), key=lambda x: x[0])
    labels_s, handles_s = zip(*pairs)
    plt.legend(handles_s, labels_s, fontsize=6, markerscale=1, bbox_to_anchor=(1.02,1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'2D_{name}_no_labels.png', dpi=1200)
    plt.close()


    #2D_PCA_centroid_labels

    plt.figure(figsize=(8,6))
    for label in sorted(species_label_list):
        sub = coords_df[coords_df['SpeciesLabel']==label]
        if not sub.empty:
            plt.scatter(sub['PC1'], sub['PC2'], color=label_color_map[label], s=50, alpha=0.7)

            cx, cy = sub['PC1'].mean(), sub['PC2'].mean()

            plt.text(cx, cy, label, fontsize=4)

    plt.xlabel(f'PC1 ({var_pct[0]:.1f}%)')
    plt.ylabel(f'PC2 ({var_pct[1]:.1f}%)')
    plt.title(f'2D PCA {name} - Centroids')
    plt.legend(labels_s, fontsize=6, markerscale=1, bbox_to_anchor=(1.02,1))
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'2D_{name}_centroids.png', dpi=1200)
    plt.close()



    #2D_PCA_strain_labels

    plt.figure(figsize=(8,6))
    for idx in coords_df.index:

        x, y = coords_df.loc[idx, ['PC1','PC2']]
        label = coords_df.loc[idx, 'SpeciesLabel']
        plt.scatter(x, y, color=label_color_map[label], s=30)
        plt.text(x+0.01, y+0.01, idx, fontsize=3)
    plt.xlabel(f'PC1 ({var_pct[0]:.1f}%)')
    plt.ylabel(f'PC2 ({var_pct[1]:.1f}%)')
    plt.title(f'2D PCA {name} - Strain Labels')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'2D_{name}_strain_labels.png', dpi=1200)
    plt.close()



    #Pairplot_PC1–PC4

    rename_map = {pcs[i]: f"{pcs[i]} ({var_pct[i]:.1f}%)" for i in range(4)}
    df_pair = coords_df.rename(columns=rename_map)
    sns.pairplot(
        df_pair, vars=list(rename_map.values()), hue='SpeciesLabel', palette=label_color_map,
        plot_kws={'alpha':0.7,'s':40}, corner=True
    )
    plt.suptitle(f'Pairwise PCs (1–4) {name}', y=1.02)
    plt.tight_layout()
    plt.savefig(f'pairplot_{name}.png', dpi=1200)
    plt.close()
    print(f"Completed {name}: PC1={var_pct[0]:.1f}%, PC2={var_pct[1]:.1f}%")



#Execute

for subset in ['Dickeya','Pectobacterium','Both']:
    species_pca(subset, subsets[subset])

print('VFDB‐based species PCA clustering complete.')

