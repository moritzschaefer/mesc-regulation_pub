import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from moritzsphd.data import ensembl_release, gn2id, targetscan_summary
from moritzsphd.util.mirna.clustering import mirbase_clusters
from scipy.stats import kstest

clusters = mirbase_clusters()

fig, ax = plt.subplots(figsize=(5.2, 3.2))
fig2, ax2 = plt.subplots(figsize=(5.2, 3.2))
fig3, ax3 = plt.subplots(figsize=(5.2, 3.2))

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

sns.ecdfplot(target_deg, ax=ax, label=f'Integrative approach (n={len(target_deg)})', color=snakemake.params['sample_colors'][mutant])
sns.ecdfplot(target_deg, ax=ax2, label=f'Integrative approach (n={len(target_deg)})', color=snakemake.params['sample_colors'][mutant])
sns.ecdfplot(target_deg, ax=ax3, label=f'Integrative approach (n={len(target_deg)})', color=snakemake.params['sample_colors'][mutant])
pd.Series({'num_targets': len(cluster_targets_minscore), 'num_up_targets': (target_deg > 0).sum(), 'ratio': (target_deg > 0).sum()/len(cluster_targets_minscore)}).to_csv(snakemake.output['up_percentage'], header=False)

# controls.extend(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist())

n_genes = len(target_deg)
assert n_genes == 324, 'miR-290 targets should be 324'

if snakemake.params['plot_ts_predictions']:
    ts_df = targetscan_summary()
    # select mir-290-295 targets with seed AAGUGCU (strongly expressed 291,294,295 members)
    ts_df = ts_df.loc[(ts_df['Species ID'] == 10090) & (ts_df['miRNA family'] == 'AAGUGCU') & ts_df['Representative miRNA'].str.match('mmu-miR-29[0-5]')].copy()
    def _get_gene_id(row):
        try:
            return ensembl_release.transcript_by_id(row['Transcript ID']).gene_id
        except ValueError:
            return gn2id(row['Gene Symbol'])
    ts_df['gene_id'] = ts_df.apply(_get_gene_id, axis=1)

    # filter the 324 genes with lowest targetscan scores (lower is "better")
    ts_target_genes = ts_df.groupby('gene_id')['Cumulative weighted context++ score'] \
                           .min().sort_values().iloc[:n_genes].index
    ts_target_genes_deg = diffexp_df.loc[ts_target_genes, 'log2FoldChange']

    sns.ecdfplot(ts_target_genes_deg, ax=ax, label=f'TargetScan predictions (n={len(ts_target_genes_deg)})', color='Gold')
    print(f'Kolmogorov-Smirnov test between DEGs of top {len(ts_target_genes_deg)} TargetScan predictions vs the predictions from the paper\'s integrative approach: {kstest(ts_target_genes_deg, target_deg)}')

if snakemake.params['plot_heap_predictions']:
    # Load HEAP peaks in 3p UTRs with miR-290-295 binding site match
    heap_df = pd.read_csv(snakemake.input['ago2_heap'], index_col=0)
    # filter for peaks in 3' UTRs and peaks with seed matches to mir-291,294,295
    heap_df = heap_df.loc[heap_df['mirna'].isin(ts_df['Representative miRNA'].drop_duplicates()) & (heap_df.is_3putr == 1)]
    heap_target_genes = heap_df.groupby('gene_id')['score'].max().sort_values(ascending=False).iloc[:n_genes].index
    heap_target_genes_deg = diffexp_df.loc[heap_target_genes, 'log2FoldChange']
    sns.ecdfplot(heap_target_genes_deg, ax=ax, label=f'AGO2-binding (n={len(heap_target_genes_deg)})', color='Brown')

    print(f'Kolmogorov-Smirnov test between DEGs of top {len(heap_target_genes_deg)} HEAP predictions vs the predictions from the paper\'s integrative approach: {kstest(heap_target_genes_deg, target_deg)}')

if snakemake.params['plot_upregulation_predictions']:
    with open(snakemake.input.up_genes, 'r') as f:
        up_genes = f.read().strip('\n').split(',')

    sns.ecdfplot(diffexp_df.loc[up_genes, "log2FoldChange"], ax=ax, label=f'RNAi_KO-upregulated (n={len(up_genes)}) ', color='#80edc5')
    print(f'Kolmogorov-Smirnov test between DEGs of {len(up_genes)} "commonly upregulated genes" vs the predictions from the paper\'s integrative approach: {kstest(diffexp_df.loc[up_genes, "log2FoldChange"], target_deg)}')

    # Supp. figure: Filtered up-genes
    up_genes_filtered = pd.Index(up_genes).intersection(ts_df.gene_id)  # This is equivalent to mir-291/4/5-seed match search in 3p-UTR
    sns.ecdfplot(diffexp_df.loc[up_genes_filtered, "log2FoldChange"], ax=ax2,
                 label=f'filtered RNAi_KO-upregs (n={len(up_genes_filtered)})', color='#60cda5')
    print(f'Kolmogorov-Smirnov test between DEGs of {len(up_genes_filtered)} "commonly upregulated genes" filtered for miR-291/294/295 binding sites (seed match AAGUGCU) vs the predictions from the paper\'s integrative approach: {kstest(diffexp_df.loc[up_genes_filtered, "log2FoldChange"], target_deg)}')

colors = snakemake.params['sample_colors'][mutant]
# Supp. figure: Interaction score has potential to improve predictions
unfiltered_interaction_df = pd.read_csv(snakemake.input['unfiltered_interaction_data']).query(f'`WT miRNA expression` > {snakemake.params.mirna_threshold}')
# Use 'Interaction score' for ranking miRNA target genes with the simple formulae: max(interaction_scores) + (log10(n) + 1) * (sum(interactions_scores)/n)
score_filtered = unfiltered_interaction_df.groupby('Geneid')['Interaction score'].agg(lambda group: group.max() + (np.log10(len(group)) + 1) * group.mean()).sort_values(ascending=False).iloc[:n_genes]
sns.ecdfplot(diffexp_df.loc[score_filtered.index, "log2FoldChange"], ax=ax2, label=f'Interaction score-based filter (n={len(score_filtered)})', color='#0000FF')
print(f'Kolmogorov-Smirnov test between DEGs of {len(score_filtered)} Interaction-score-filtered targets vs the predictions from the paper\'s integrative approach: {kstest(diffexp_df.loc[score_filtered.index, "log2FoldChange"], target_deg)}')

# plot low_up_genes
with open(snakemake.input.low_up_genes, 'r') as f:
    low_up_genes = f.read().strip('\n').split(',')
sns.ecdfplot(diffexp_df.loc[low_up_genes, "log2FoldChange"], ax=ax3, label=f'lowly upregulated genes (n={len(low_up_genes)})', color='#7ec0ee')

low_predictions = pd.read_csv(snakemake.input['low_mirna_targets'], index_col=[0, 1, 2, 3, 4, 5])
low_cluster_predictions = low_predictions.loc[low_predictions.apply(lambda row: row.name[5] in target_cluster, axis=1)]
low_cluster_targets = low_cluster_predictions.groupby(low_cluster_predictions.index.get_level_values(0))['Interaction score'].sum().sort_values(ascending=False)  # irrelevant
low_target_deg = diffexp_df.loc[low_cluster_targets.index, 'log2FoldChange']
sns.ecdfplot(low_target_deg, ax=ax3, label=f'Low-up genes + Integrative approach (n={len(low_target_deg)})', color='#Ff7f50')

# just take the last one as control and format plot a little
for axis in [ax, ax2, ax3]:
    sns.ecdfplot(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist(), ax=axis, label='expressed genes (control)', color='gray')

    axis.set_xlim([-0.5, 1.6])
    axis.legend(loc='lower right')
    axis.axhline(0.5, color='gray')
    axis.axvline(0, color='gray')
    axis.grid(False)
    axis.set_xlabel('log2FoldChange in miR-290-295_KO vs WT')
    axis.set_title('Differential expression of predicted miR-290-295 targets')
    sns.despine(ax=axis)
# plt.legend(loc='lower left', bbox_to_anchor=(0.1, 0.05))
# plt.legend(loc='lower right')
fig.savefig(snakemake.output['plot'])
fig2.savefig(snakemake.output['supp_plot'])
fig3.savefig(snakemake.output['low_up_supp_plot'])
