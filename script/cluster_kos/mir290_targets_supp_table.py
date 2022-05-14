import pandas as pd
from moritzsphd.util.mirna.clustering import mirbase_clusters

df = pd.read_csv(snakemake.input['interactions'])
cluster_ko_des = pd.read_excel(snakemake.input['mir_cluster_des'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1]).xs('miR-290-295', axis=1, level=0).droplevel(level=1)

clusters = mirbase_clusters()
mir295mirnas = [c for c in clusters if 'mmu-miR-295-5p' in c][0]

mir290_targets = df.loc[df['miRNA'].isin(mir295mirnas), ['Geneid', 'Gene name']].set_index('Geneid').drop_duplicates()

mir290_targets['miR-290-295_KO log2FC'] = cluster_ko_des['log2FoldChange']
mir290_targets['miR-290-295_KO padj'] = cluster_ko_des['padj']
mir290_targets['is TF'] = 'no'
tf_df = pd.read_csv(snakemake.input['tf_annotation'], index_col=0)
mir290_targets.loc[mir290_targets.index.intersection(tf_df.index), 'is TF'] = 'yes'

mir290_targets['up/down'] = mir290_targets['miR-290-295_KO log2FC'].apply(lambda v: 'up' if v > 0 else 'down')
mir290_targets['Statistically significant'] = mir290_targets.apply(lambda v: 'yes' if v['miR-290-295_KO log2FC'] > 0 and v['miR-290-295_KO padj'] < 0.1 else 'no', axis=1)

mir290_targets['Status'] = 'Unconfirmed'
mir290_targets.loc[(mir290_targets['miR-290-295_KO log2FC'] > 0.0),
                   'Status'] = 'Positive log2FC'
mir290_targets.loc[(mir290_targets['miR-290-295_KO log2FC'] > snakemake.params['log2fc_threshold']) &
                   (mir290_targets['miR-290-295_KO padj'] < snakemake.params['padj_threshold']),
                   'Status'] = 'Significantly Upregulated/Confirmed'

mir290_targets.to_excel(snakemake.output['supp_table'])
