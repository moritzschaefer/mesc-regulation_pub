import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from moritzsphd.util.mirna.clustering import mirbase_clusters

df = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)
# df = np.log2(df+1)

clusters = mirbase_clusters()
mutant = 'miR-290-295'
mirna_id = mutant[3:mutant.find('_')] if '_' in mutant else mutant[4:7]
target_cluster = [c for c in clusters if f'mmu-miR-{mirna_id}-3p' in c or f'mmu-miR-{mirna_id}a-3p' in c][0]

fig, ax = plt.subplots(figsize=(3.5, 3.5))

# plt.plot([0, df[['WT_1', 'WT_2']].max().max()], [0, df[['WT_1', 'WT_2']].max().max()], color='gray')
df['mirna_group'] = 'other'
df.loc[df.index.intersection(target_cluster), 'mirna_group'] = 'miR-290-295'
ax.plot([1, df[['RIP_AGO2', 'RIP_AGO1']].max().max()], [1, df[['RIP_AGO2', 'RIP_AGO1']].max().max()], color='gray')
ax = sns.scatterplot(data=df.sort_values('mirna_group', ascending=False), x='RIP_AGO2', y='RIP_AGO1', hue='mirna_group', palette='pastel', s=50)
txts = []
for mirna in target_cluster:
    try:
        txts.append(ax.text(df.loc[mirna, 'RIP_AGO2'], df.loc[mirna, 'RIP_AGO1'], mirna))
    except KeyError:
        pass

ax.set_xscale('log')
ax.set_yscale('log')
sns.despine()
adjust_text(txts, arrowprops=dict(arrowstyle="->", color='black', lw=0.5))
ax.set_yticks([10, 1000, 100000])
ax.set_xticks([10, 1000, 100000])

plt.savefig(snakemake.output[0])
