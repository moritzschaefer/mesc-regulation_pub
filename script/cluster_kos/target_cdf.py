import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from moritzsphd.data import ensembl_release, gn2id, targetscan_summary
from moritzsphd.util.mirna.clustering import mirbase_clusters
from scipy.stats import kstest

clusters = mirbase_clusters()

fig, ax = plt.subplots(figsize=(3.0, 2.0))
fig2, ax2 = plt.subplots(figsize=(3.0, 2.0))

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

# alternative control
non_targets = set(predictions.index.get_level_values(0)) - set(cluster_targets_minscore)
non_target_deg = diffexp_df.loc[non_targets, 'log2FoldChange']

sns.ecdfplot(target_deg, ax=ax, label=f'Integrative approach (n={len(target_deg)})', color=snakemake.params['sample_colors'][mutant])
sns.ecdfplot(target_deg, ax=ax2, label=f'Integrative approach (n={len(target_deg)})', color=snakemake.params['sample_colors'][mutant])
pd.Series({'num_targets': len(cluster_targets_minscore), 'num_up_targets': (target_deg > 0).sum(), 'ratio': (target_deg > 0).sum()/len(cluster_targets_minscore)}).to_csv(snakemake.output['up_percentage'], header=False)
print(f'Mean target DEG: {target_deg.mean()}')

# controls.extend(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist())

n_genes = len(target_deg)

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

    # filter the <n_genes> genes with lowest targetscan scores (lower is "better")
    ts_target_genes = ts_df.groupby('gene_id')['Cumulative weighted context++ score'] \
                           .min().sort_values().iloc[:n_genes].index
    ts_target_genes_deg = diffexp_df.loc[ts_target_genes, 'log2FoldChange']

    sns.ecdfplot(ts_target_genes_deg, ax=ax, label=f'TargetScan predictions (n={len(ts_target_genes_deg)})', color='Gold')
    print(f'Mean TargetScan DEG: {ts_target_genes_deg.mean()}')
    print(f'Kolmogorov-Smirnov test between DEGs of top {len(ts_target_genes_deg)} TargetScan predictions vs the predictions from the paper\'s integrative approach: {kstest(ts_target_genes_deg, target_deg)}')

if snakemake.params['plot_heap_predictions']:
    # Load HEAP peaks in 3p UTRs with miR-290-295 binding site match
    heap_df = pd.read_csv(snakemake.input['ago2_heap'], index_col=0)
    # filter for peaks in 3' UTRs and peaks with seed matches to mir-291,294,295
    heap_df = heap_df.loc[heap_df['mirna'].isin(ts_df['Representative miRNA'].drop_duplicates()) & (heap_df.is_3putr == 1)]
    heap_target_genes = heap_df.groupby('gene_id')['score'].max().sort_values(ascending=False).iloc[:n_genes].index
    heap_target_genes_deg = diffexp_df.loc[heap_target_genes, 'log2FoldChange']
    sns.ecdfplot(heap_target_genes_deg, ax=ax, label=f'AGO2-binding (n={len(heap_target_genes_deg)})', color='Brown')

    print(f'Mean HEAP DEG: {heap_target_genes_deg.mean()}')
    print(f'Kolmogorov-Smirnov test between DEGs of top {len(heap_target_genes_deg)} HEAP predictions vs the predictions from the paper\'s integrative approach: {kstest(heap_target_genes_deg, target_deg)}')

if snakemake.params['plot_upregulation_predictions']:
    with open(snakemake.input.up_genes, 'r') as f:
        up_genes = f.read().strip('\n').split(',')
    up_genes_deg = diffexp_df.loc[up_genes, 'log2FoldChange']

    sns.ecdfplot(up_genes_deg, ax=ax, label=f'RNAi_KO-upregulated (n={len(up_genes)}) ', color='#80edc5')
    print(f'Mean upreg DEG: {up_genes_deg.mean()}')
    print(f'Kolmogorov-Smirnov test between DEGs of {len(up_genes)} "commonly upregulated genes" vs the predictions from the paper\'s integrative approach: {kstest(up_genes_deg, target_deg)}')

colors = snakemake.params['sample_colors'][mutant]

# plot low_up_genes
with open(snakemake.input.low_up_genes, 'r') as f:
    low_up_genes = f.read().strip('\n').split(',')
low_up_genes_deg = diffexp_df.loc[low_up_genes, "log2FoldChange"]
print(f'Mean for low_up_genes_deg: {low_up_genes_deg.mean()}')
sns.ecdfplot(low_up_genes_deg, ax=ax2, label=f'lowly upregulated genes (n={len(low_up_genes)})', color='#7ec0ee')

low_predictions = pd.read_csv(snakemake.input['low_mirna_targets'], index_col=[0, 1, 2, 3, 4, 5])
low_cluster_predictions = low_predictions.loc[low_predictions.apply(lambda row: row.name[5] in target_cluster, axis=1)]
low_cluster_targets = low_cluster_predictions.groupby(low_cluster_predictions.index.get_level_values(0))['Interaction score'].sum().sort_values(ascending=False)  # irrelevant
low_target_deg = diffexp_df.loc[low_cluster_targets.index, 'log2FoldChange']
print(f'Mean for integrative-filtered low_up_genes_deg: {low_target_deg.mean()}')
sns.ecdfplot(low_target_deg, ax=ax2, label=f'Low-up genes + Integrative approach (n={len(low_target_deg)})', color='#Ff7f50')

# other control (can be shown optionally)
# sns.ecdfplot(non_target_deg, ax=ax, label=f'non-miR-290-295 predicted target, (n={len(non_target_deg)})', color='black')

# just take the last one as control and format plot a little
for f, axis, output_label in zip([fig, fig2], [ax, ax2], ['plot', 'low_up_supp_plot']):
    sns.ecdfplot(diffexp_df.loc[diffexp_df.baseMean > 100, 'log2FoldChange'].tolist(), ax=axis, label='expressed genes (control)', color='gray')

    axis.set_xlim([-0.5, 1.6])
    axis.legend(loc='lower right')
    axis.axhline(0.5, color='gray')
    axis.axvline(0, color='gray')
    axis.grid(False)
    axis.set_xlabel('log2FoldChange in miR-290-295_KO vs WT')
    axis.set_title('Diff. expression of predicted miR-290-295 targets')
    sns.despine(ax=axis)
    f.set_tight_layout(True)

    f.savefig(snakemake.output[output_label])
# plt.legend(loc='lower left', bbox_to_anchor=(0.1, 0.05))
# plt.legend(loc='lower right')
