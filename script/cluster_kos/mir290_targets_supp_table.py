import pandas as pd
from moritzsphd.util.mirna.clustering import mirbase_clusters

df = pd.read_excel(snakemake.input['interactions'], sheet_name='Filtering 4', skiprows=2, index_col=[0, 1, 2])
cluster_ko_des = pd.read_excel(snakemake.input['mir_cluster_des'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1]).xs('miR-290-295', axis=1, level=0).droplevel(level=1)

clusters = mirbase_clusters()
mir295mirnas = [c for c in clusters if 'mmu-miR-295-5p' in c][0]

mir290_targets = df.loc[df.index.get_level_values('miRNA').isin(mir295mirnas)]

mir290_targets = mir290_targets.index.to_frame().set_index('Geneid')[['Gene name']].drop_duplicates()

mir290_targets['miR-290-295_KO log2FC'] = cluster_ko_des['log2FoldChange']
mir290_targets['miR-290-295_KO padj'] = cluster_ko_des['padj']
mir290_targets['is TF'] = 'no'
tf_df = pd.read_csv(snakemake.input['tf_annotation'], index_col=0)
mir290_targets.loc[mir290_targets.index.intersection(tf_df.index), 'is TF'] = 'yes'


mir290_targets['Status'] = 'Unconfirmed'
mir290_targets.loc[(mir290_targets['miR-290-295_KO log2FC'] > 0.0),
                   'Status'] = 'Positive log2FC'
mir290_targets.loc[(mir290_targets['miR-290-295_KO log2FC'] > 0.5) &
                   (mir290_targets['miR-290-295_KO padj'] < 0.1),
                   'Status'] = 'Significantly Upregulated/Confirmed'

mir290_targets.to_excel(snakemake.output['supp_table'])
