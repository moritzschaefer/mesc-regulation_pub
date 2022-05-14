import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from moritzsphd.data import (ensembl_release, gn2id, grc_annotation,
                             remote_file, targetscan_summary)
from moritzsphd.util.mirna.clustering import mirbase_clusters
from scipy.stats import gaussian_kde, kstest, pearsonr

######## HEAP MRES ##############
# take unfiltered interactions
# from moritzsphd.data import mrna_expression, average_sample, gid2n

df = pd.read_csv(snakemake.input.unfiltered_interactions)
gene_max_scores = df.groupby('Geneid')['Interaction score'].max()

# filter by min_mirna_expression and TargetScan score > 0.0
filtered = df[(df['WT miRNA loading'] > snakemake.params['mirna_threshold']) & (df['weighted context++ score'] < 0)]
gene_interaction_counts = filtered.groupby('Geneid')['miRNA'].count()
gene_interaction_counts.name = 'num_ts_interactions'
gene_interaction_counts = gene_interaction_counts.to_frame()

# also take ago2_heap and get number of peaks and mean peak size from heap_ago2 on a per_gene basis. Now correlate on a per-gene basis the number of interactions and the peak values
heap_df = pd.read_csv(snakemake.input.ago2_heap_peaks, index_col=None)
def _find_target(row):
    hits = grc_annotation().region(seqid=row['seqnames'].replace('chr', ''),
                                   start=row['start'], end=row['end'], strand=row['strand'])
    return [hit.attributes['gene_id'][0] for hit in hits]
# hits = pd.DataFrame([{
#     'gene_id': hit.attributes['gene_id'][0],
#     'protein_coding': 'protein_coding' in hit.attributes['gene_biotype'],
#     'feature_type': hit.featuretype
#     } for hit in hits])

gene_ids = heap_df.apply(_find_target, axis=1)
heap_df['gene_id'] = gene_ids.apply(np.unique)
heap_df = heap_df.explode('gene_id')
# show mean heap signal in function of number of interactions
# show ribo DEG also for groups
# gene_ids = gene_interaction_counts.index.intersection(num_heap_peaks.index)
gene_interaction_counts['num_heap_peaks'] = heap_df.groupby('gene_id')['score'].count()
gene_interaction_counts['mean_heap_signal_score'] = heap_df.groupby('gene_id')['score'].mean()

interactions = pd.read_csv(snakemake.input.interaction_ranking)

fig, ax = plt.subplots(figsize=(1.5, 2.7))

ts_bin_thres = snakemake.params['ts_mre_count_bins'] + [gene_interaction_counts['num_ts_interactions'].max()+1]
ts_gene_groups = {f'{b1}-{b2-1}': gene_interaction_counts.index[(gene_interaction_counts.num_ts_interactions >= b1) & (gene_interaction_counts.num_ts_interactions < b2)]
                  for b1, b2 in zip(ts_bin_thres[:-1], ts_bin_thres[1:])}

# group
gene_interaction_counts['interaction_count'] = gene_interaction_counts.index.map(lambda gid: [i for i, g in ts_gene_groups.items() if gid in g][0])

sns.boxplot(data=gene_interaction_counts, x='interaction_count', y='num_heap_peaks', order=ts_gene_groups.keys(), ax=ax, color='gray')
ax.set_ylim([0, 15])
_ = ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')
fig.savefig(snakemake.output.heap_for_num_mres)


##################### by ribo-seq ##################

interaction_count = interactions.groupby('Geneid').apply(len)

ribo_df = pd.read_excel(snakemake.input['ribo_data'], index_col=[0, 1], skiprows=2)
ribo_df.rename(columns={'AGO2&1_log2FC': 'AGO12_log2FC', 'AGO2&1_padj': 'AGO12_padj'}, inplace=True)
mutants = snakemake.params['mutants']
ribo_log2fc = ribo_df[[m.upper() + '_log2FC' for m in mutants]]
ribo_log2fc.columns = mutants
ribo_log2fc.reset_index(level=1, inplace=True, drop=True)
# show for all three together
interaction_count_ribo = interaction_count.loc[ribo_log2fc.index.intersection(interaction_count.index)]

number_gene_interactions_count = interactions.groupby('Geneid').count()['Gene name']
bin_thres = snakemake.params['mre_count_bins'] + [max(number_gene_interactions_count) + 1]
ribo_gene_groups = {f'{b1}-{b2-1}': interaction_count_ribo.index[(interaction_count_ribo >= b1) & (interaction_count_ribo < b2)]
               for b1, b2 in zip(bin_thres[:-1], bin_thres[1:])}
ribo_misreg = ribo_log2fc.loc[interaction_count_ribo.index]
ribo_misreg['interaction_count'] = ribo_misreg.index.map(lambda gid: [i for i, g in ribo_gene_groups.items() if gid in g][0])
fig, ax = plt.subplots(figsize=(3.5, 2.7))

for mutant in mutants:
    # color = snakemake.params.sample_colors[mutant]
    # sns.regplot(x=ribo_misreg[mutant], y=interaction_count_ribo, scatter_kws={'s': 5, 'alpha': 0.5}, color=color, ax=ax, label=f'{mutant}_KO')
    pearsonr(ribo_misreg[mutant], interaction_count_ribo)
    # sns.despine()
    print(f'{mutant} pearsonr: {pearsonr(ribo_misreg[mutant], interaction_count_ribo)[0]:.2}')

sns.boxplot(data=ribo_misreg.melt(id_vars=['interaction_count'], value_vars=mutants, value_name='log2FC', var_name='mutant'), x='mutant', y='log2FC', hue='interaction_count', ax=ax, hue_order=ribo_gene_groups.keys())
# ax.set_title('RNAi_KO Ribosome occupancy in function of number of targeting miRNAs')
sns.despine()
ax.set_xlabel('Diff. ribosome occupancy (log2FC)')
ax.set_ylabel('Number of miRNA interations')
ax.set_title('Number of targeting miRNAs vs. RNAi_KO ribo occupancy for each target')
ax.legend()
ax.set_ylim([ax.get_ylim()[0], 4.0])
fig.savefig(snakemake.output.interaction_count_ribo)


############### mir-290-signal ##############

clusters = mirbase_clusters()
mutant = 'miR-290-295'
mirna_id = mutant[3:mutant.find('_')] if '_' in mutant else mutant[4:7]
target_cluster = [c for c in clusters if f'mmu-miR-{mirna_id}-3p' in c or f'mmu-miR-{mirna_id}a-3p' in c][0]

diffexp_df = pd.read_excel(snakemake.input['quantseq_data'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1])
diffexp_df = diffexp_df.xs(mutant, axis=1, level=0).droplevel(level=1)

cluster_predictions = interactions.loc[interactions['miRNA'].isin(target_cluster)]

interaction_count_mir290 = cluster_predictions.groupby('Geneid').apply(len)

fig, ax = plt.subplots(figsize=(1.5, 2.3))

mir290_gene_groups = {f'{b1}-{b2-1}': interaction_count_mir290.index[(interaction_count_mir290 >= b1) & (interaction_count_mir290 < b2)] for b1, b2 in zip(bin_thres[:-1], bin_thres[1:])}

mir290_misreg = diffexp_df.loc[interaction_count_mir290.index]
mir290_misreg['interaction_count'] = mir290_misreg.index.map(lambda gid: [i for i, g in mir290_gene_groups.items() if gid in g][0])

sns.boxplot(data=mir290_misreg, x='interaction_count', y='log2FoldChange', order=mir290_gene_groups.keys(), ax=ax)
sns.despine()
ax.set_title('miR-290 targets in miR-290_KO')
_ = ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha='right')

ax.set_ylim([ax.get_ylim()[0], 2.0])
plt.tight_layout()
fig.savefig(snakemake.output.interaction_count_mir290)
