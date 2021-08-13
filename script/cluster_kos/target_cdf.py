import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from moritzsphd.data import ensembl_release, gn2id, targetscan_summary
from moritzsphd.util.mirna.clustering import mirbase_clusters

clusters = mirbase_clusters()
fig, ax = plt.subplots(figsize=(4.7, 3.2))

# controls = []

predictions = pd.read_csv(snakemake.input['mirna_targets'], index_col=[0, 1, 2, 3, 4, 5])
mutant = snakemake.wildcards['cluster']

diffexp_df = pd.read_excel(snakemake.input['data'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1])
diffexp_df = diffexp_df.xs(mutant, axis=1, level=0).droplevel(level=1)

mirna_id = mutant[3:mutant.find('_')] if '_' in mutant else mutant[4:7]
target_cluster = [c for c in clusters if f'mmu-miR-{mirna_id}-3p' in c or f'mmu-miR-{mirna_id}a-3p' in c][0]
cluster_predictions = predictions.loc[predictions.apply(lambda row: row.name[5] in target_cluster, axis=1)]

cluster_targets = cluster_predictions.groupby(cluster_predictions.index.get_level_values(0))['Interaction score'].sum().sort_values(ascending=False)

# cluster_targets_minscore = cluster_targets.index[cluster_targets > 2]
cluster_targets_minscore = cluster_targets.index
target_deg = diffexp_df.loc[cluster_targets_minscore, 'log2FoldChange']

sns.ecdfplot(target_deg, ax=ax, label=f'Predicted miR-290-295 targets', color=snakemake.params['sample_colors'][mutant])
pd.Series({'num_targets': len(cluster_targets_minscore), 'num_up_targets': (target_deg > 0).sum(), 'ratio': (target_deg > 0).sum()/len(cluster_targets_minscore)}).to_csv(snakemake.output['up_percentage'], header=False)

# controls.extend(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist())


if snakemake.params['plot_ts_predictions']:
    ts_df = targetscan_summary()
    ts_df = ts_df.loc[(ts_df['Species ID'] == 10090) & ts_df['Representative miRNA'].str.match('mmu-miR-29[0-5]') & (ts_df['Cumulative weighted context++ score'] < -0.4)].copy()

    def _get_gene_id(row):
        try:
            return ensembl_release.transcript_by_id(row['Transcript ID']).gene_id
        except ValueError:
            return gn2id(row['Gene Symbol'])
    ts_df['gene_id'] = ts_df.apply(_get_gene_id, axis=1)
    sns.ecdfplot(diffexp_df.loc[ts_df['gene_id'].dropna().unique(), 'log2FoldChange'], ax=ax, label='TargetScan miR-290-295 predictions', color='Gold')

if snakemake.params['plot_heap_predictions']:
    heap_df = pd.read_csv(snakemake.input['ago2_heap'], index_col=0)
    # get top 10000 peaks/interactions and extract genes
    genes = heap_df.sort_values('score', ascending=False).iloc[:10000]['gene_id'].unique()
    sns.ecdfplot(diffexp_df.loc[genes, 'log2FoldChange'], ax=ax, label='miR-290-295 AGO2-binding targets', color='Brown')

if snakemake.params['plot_upregulation_predictions']:
    with open(snakemake.input.up_genes, 'r') as f:
        up_genes = f.read().strip('\n').split(',')
    sns.ecdfplot(diffexp_df.loc[up_genes, 'log2FoldChange'], ax=ax, label='RNAi mutants-upregulated genes', color='#80edc5')

# just take the last one as control
sns.ecdfplot(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist(), ax=ax, label='expressed genes (control)', color='gray')

ax.set_xlim([-1, 2])
ax.legend(loc='lower right')
ax.axhline(0.5, color='gray')
ax.axvline(0, color='gray')
ax.grid(False)
ax.set_title('Differential expression in miR-290-295 KO')
sns.despine()
# plt.legend(loc='lower left', bbox_to_anchor=(0.1, 0.05))
plt.legend(loc='upper left')
fig.savefig(snakemake.output['plot'])
