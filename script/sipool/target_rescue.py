import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from moritzsphd.data import gid2n, gn2id  # we anyways just need tfap4..
from moritzsphd.util.mirna.clustering import mirbase_clusters

mir290de = pd.read_excel(snakemake.input['cluster_ko'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1]).xs('miR-290-295', axis=1, level=0).droplevel(level=1, axis=0)
clusters = mirbase_clusters()
mir295mirnas = [c for c in clusters if 'mmu-miR-295-5p' in c][0]

interactions = pd.read_csv(snakemake.input['interactions'])
mir290_targets = interactions.loc[interactions.miRNA.isin(mir295mirnas), 'Geneid'].drop_duplicates()
mir290de_genes = mir290de.loc[(mir290de.log2FoldChange.abs() > snakemake.params['log2fc_threshold']) & (mir290de.padj < snakemake.params['padj_threshold'])]
mir290de_genes['role'] = pd.Categorical(['unexplained'] * len(mir290de_genes), categories=['unexplained', 'miR290 target', 'potential TF target'], ordered=True)
mir290de_genes.loc[mir290de_genes.index.isin(mir290_targets) & (mir290de_genes.log2FoldChange > 0), 'role'] = 'miR290 target'
genes = snakemake.params['genes']
out = []
acceptance_factor = snakemake.params['acceptance_factor']
fig, axes = plt.subplots(1, len(genes), figsize=(4 * len(genes), 4), sharey=True, squeeze=False)

mir290de_genes.rename(columns={'log2FoldChange': 'miR290KO log2FC'}, inplace=True)
#skiprows=2 is relevant if the xlsx got beautified (for publication)
si_des = pd.read_excel(snakemake.input['sipool_des'], sheet_name='Main', skiprows=1,
                       index_col=[0, 1], header=[0, 1]).droplevel(level=1, axis=0)

for gene, ax in zip(genes, axes[0, :]):
    sidf = si_des.xs(f'si_{gene.lower()}', axis=1, level=0)
    plotdf = mir290de_genes.copy()
    plotdf['siPOOL log2FC'] = sidf['log2FoldChange']

    ratio = plotdf['miR290KO log2FC'] / plotdf['siPOOL log2FC']
    selector = (ratio > -acceptance_factor) & (ratio < -1/acceptance_factor)
    mir290de_genes[f'si{gene}_ratio'] = ratio
    mir290de_genes[f'si{gene}_rescued'] = False
    mir290de_genes.loc[selector, f'si{gene}_rescued'] = True

    plotdf.loc[selector, 'role'] = 'potential TF target'

    ax.axline((0, 0), slope=(-acceptance_factor), color='gray')
    ax.axline((0, 0), slope=(-1/acceptance_factor), color='gray')

    sns.scatterplot(data=plotdf.sort_values('role'),
                    x='miR290KO log2FC', y='siPOOL log2FC', hue='role',
                    palette={'potential TF target': 'orange', 'unexplained': 'gray', 'miR290 target': '#6ca6cd'},
                    s=15, linewidth=0, edgecolor='none',
                    rasterized=True,
                    legend=(gene == genes[0]),
                    ax=ax
                    )
    ax.plot(plotdf.loc[gn2id(gene), 'miR290KO log2FC'], plotdf.loc[gn2id(gene), 'siPOOL log2FC'],'ro')
    ax.axhline(0, color='black')
    ax.axvline(0, color='black')
    ax.set_title(f'siPOOL for {gene}')
    ax.set_xlim([-5, 5])
    ax.set_ylim([-5, 5])
    ax.text(-4, 3, str((selector & (plotdf['siPOOL log2FC'] > 0)).sum()))
    ax.text(3, -3, str((selector & (plotdf['siPOOL log2FC'] < 0)).sum()))
    sns.despine()

fig.tight_layout()
fig.savefig(snakemake.output['plot_png'])
fig.savefig(snakemake.output['plot_svg'])
mir290de_genes['Gene name'] = mir290de_genes.index.map(gid2n)
mir290de_genes.index.name = 'Geneid'
mir290de_genes.set_index('Gene name', append=True, inplace=True)
mir290de_genes.to_csv(snakemake.output['rescue_df'])
