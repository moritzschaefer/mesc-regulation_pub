import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from moritzsphd.data import gn2id

sidf = pd.read_csv(snakemake.input['quant_seq'], index_col=0, sep='\t')
sidf *= 1e6/sidf.sum()  # CPM normalization

mir290de = pd.read_excel(snakemake.input['cluster_ko'], sheet_name='Main', skiprows=,, index_col=[0, 1], header=[0, 1]).xs('miR-290-295', axis=1, level=0).droplevel(level=1, axis=0)
mir290de = mir290de.loc[[gn2id(g) for g in snakemake.params['genes']], ['log2FoldChange']]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3.5, 1.7), gridspec_kw={'width_ratios': [1, 5]}, sharey=True)
cbar_ax = fig.add_axes([.91, .3, .03, .4])

# ax1 =  # plot here somehow the upregulation in miR290

# normalize to si_neg per replicate
for column in sidf.columns:
    if column.startswith('si_neg'):
        continue
    sidf[column] /= sidf[f'si_neg_{column[-1]}']

sidf.drop(columns=['si_neg_1', 'si_neg_2'], inplace=True)

# log2
sidf = np.log2(sidf)
sidf = sidf.loc[[gn2id(g) for g in snakemake.params['genes']]]

vmin = min(sidf.min().min(), mir290de.min().min())
vmax = max(sidf.max().max(), mir290de.max().max())

sns.heatmap(mir290de,
            ax=ax1,
            cbar=True,
            annot=True,
            cmap='coolwarm',
            center=0,
            vmin=vmin, vmax=vmax,
            cbar_kws={'label': 'log2-foldchange'},
            cbar_ax=cbar_ax)
sns.heatmap(sidf,
            ax=ax2,
            cmap='coolwarm',
            annot=True,
            center=0,
            vmin=vmin, vmax=vmax,
            cbar=False,
            cbar_ax=None)

# cosmetics
ax2.set_yticklabels(snakemake.params['genes'])
ax2.set_xticklabels(['si' + s.get_text()[3:].title() for s in ax2.get_xticklabels()])
ax1.set_xticklabels(['miR290-295 KO vs. WT'], rotation=25, ha='right')
ax2.set_xlabel('TF-siPOOL vs neg. contr. in KO cells')
ax2.axvline(2, color='black')
ax2.axvline(4, color='black')

fig.tight_layout(rect=[0, 0, .85, 1])
fig.savefig(snakemake.output[0])
