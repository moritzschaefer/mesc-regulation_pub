import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.util.mirna.clustering import mirbase_clusters

clusters = mirbase_clusters()
mutant = 'miR-290-295'
mirna_id = mutant[3:mutant.find('_')] if '_' in mutant else mutant[4:7]
target_cluster = [c for c in clusters if f'mmu-miR-{mirna_id}-3p' in c or f'mmu-miR-{mirna_id}a-3p' in c][0]

df = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)

degs = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0, sheet_name='Differential loading')

fig, ax = plt.subplots(figsize=(6, 6))

df['differential'] = 'non-significant'
df.loc[df.index.intersection(degs.index[degs.padj < snakemake.params['padj_threshold']]), 'differential'] = 'significant'
df['Agos_RIP_mean'] = (df['RIP_AGO2'] + df['RIP_AGO1']) / 2
sns.scatterplot(data=df.sort_values('differential', ascending=False), x='Agos_RIP_mean', hue='differential', y='Input', palette='pastel', s=50, color='black', ax=ax)

plt.xscale('log')
plt.yscale('log')
sns.despine()
plt.yticks([10, 1000, 100000])
plt.xticks([10, 1000, 100000])

plt.savefig(snakemake.output[0])
