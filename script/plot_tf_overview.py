import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.util.mirna.clustering import mirbase_clusters

clusters = mirbase_clusters()
mir295mirnas = [c for c in clusters if 'mmu-miR-295-5p' in c][0]
mir290de = pd.read_excel(snakemake.input['supp_table'])
df = pd.read_csv(snakemake.input['interaction_ranking'])
df = df[df['miRNA'].isin(mir295mirnas) & df['is_tf']]
positives = mir290de.loc[(mir290de['miR-290-295_KO log2FC'] > snakemake.params['log2fc_threshold']) & (mir290de['miR-290-295_KO padj'] < snakemake.params['padj_threshold'])]

df = df[df.Geneid.isin(positives['Geneid'])]

metrics = df.groupby('Gene name')['Interaction score'].agg(['max', 'size', 'mean'])
metrics['miR-290-295 log2FC'] = positives.set_index('Gene name')['miR-290-295_KO log2FC']
x = 'Number of targeting miR-290-295-miRNAs'
y = 'Maximum interaction score'
hue = 'score mean'
metrics.rename(columns={'size': x, 'max': y, 'mean': hue}, inplace=True)
fig, ax = plt.subplots(figsize=(6, 4))

ax = sns.scatterplot(data=metrics, y=y, x=x, size=hue, hue=hue, ax=ax)
ax.set_xlabel(x)
ax.set_ylabel(y)
ax.set_xlim([0.5, 15])
sns.despine()
for index, row in metrics.iterrows():
    ax.text(row[x] + 0.15, row[y] + 0.01, index)
fig.savefig(snakemake.output[0])
