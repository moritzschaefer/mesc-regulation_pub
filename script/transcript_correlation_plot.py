import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

df = pd.read_csv(snakemake.input['expression'], sep='\t', index_col=0)
# CPM normalization
df *= (1e6/df.sum())
df = df.groupby(df.columns.map(lambda v: v[:v.find('_')]), axis=1).mean()
df = df[['WT'] + snakemake.params['mutants']]

# cut away genes that are lowly expressed in all samples
df = df[(df > 1).any(axis=1)]

# df = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
# # cut away genes that are lowly expressed in all samples
# df = df.loc[(df.xs('tpm_expression', axis=1, level=1) > 1).any(axis=1)]
# # get log2FoldChange
# df = df.xs('tpm_expression', axis=1, level=1)
# df = df[['WT'] + snakemake.params['mutants']]

fig, ax = plt.subplots(figsize=(5, 5))
cmap = sns.color_palette("flare", as_cmap=True)

sns.set_theme()

samples = [s[:s.find('_')] for s in df.columns]
# Create a categorical palette to identify the networks
sample_pal = sns.husl_palette(df.shape[1] // 2, s=.45)
sample_lut = dict(zip(map(str, np.unique(samples)), sample_pal))

# Convert the palette to vectors that will be drawn on the side of the matrix
sample_colors = pd.Series(samples, index=df.columns).map(sample_lut)

# Draw the full plot
# g = sns.clustermap(df.corr(), cmap=cmap,
#                    # row_colors=sample_colors, col_colors=sample_colors,
#                    # dendrogram_ratio=(.1, .2),
#                    # cbar_pos=(.02, .32, .03, .2),
#                    linewidths=.75, figsize=(6, 6))
# g.ax_row_dendrogram.remove()

sns.heatmap(df.corr(), cmap=cmap, square=True, linewidths=.5, cbar_kws={"shrink": .5}, ax=ax)
# sns.heatmap(df[data_cols].corr(), annot=True, cmap='coolwarm', linewidths=3, linecolor='black')

plt.savefig(snakemake.output[0])
