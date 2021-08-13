import re

import pandas as pd

from moritzsphd.data import gid2n
from moritzsphd.util.mirna.clustering import mirbase_clusters

with open(snakemake.input['model_log'], 'r') as f:
    model_log = f.read()
mean_x = float(re.search('Mean of positive log2FoldChanges: ([0-9\.]+)', model_log).groups()[0])

ml_prediction = pd.read_csv(snakemake.input['ml_prediction'])
manual_prediction = pd.read_csv(snakemake.input['manual_prediction'])
manual_prediction_grouped = pd.read_csv(snakemake.input['manual_prediction_grouped'], index_col=0)

clusters = mirbase_clusters()
mir290cluster = [c for c in clusters if 'mmu-miR-291a-3p' in c][0]


ml_prediction = ml_prediction[ml_prediction['pred_log2fc'] > mean_x].groupby('Geneid').agg({'pred_log2fc': 'max', 'miRNAs': ','.join})
ml_prediction.rename(columns={'miRNAs': 'ML_miRNAs', 'pred_log2fc': 'ML_max_pred_log2FC'}, inplace=True)

with open(snakemake.output['num_mir290_up'], 'w') as f:
    mir290_targets = ml_prediction[ml_prediction['ML_miRNAs'].apply(lambda mirnas: len(mir290cluster & set(mirnas.split(','))) > 0)]
    # now check, which of these are actually up

    mir290de = pd.read_excel(snakemake.input['mir_cluster_des'], sheet_name='Main', skiprows=2, index_col=[0, 1], header=[0, 1]).xs('miR-290-295', axis=1, level=0).droplevel(axis=0, level=1)
    mir290_up_targets = mir290de.query('log2FoldChange > 0').index.intersection(mir290_targets.index)

    f.write(str(len(mir290_up_targets)))

manual_prediction = manual_prediction.groupby('Geneid')['miRNA'].apply(lambda v: ','.join(v)).to_frame()
manual_prediction['manual_integration_score'] = manual_prediction_grouped['score']
manual_prediction.rename(columns={'miRNA': 'manual_integration_miRNAs'}, inplace=True)
df = manual_prediction.join(ml_prediction, how='outer')

df['identification_source'] = df.apply(lambda row: 'both' if (not pd.isna(row['manual_integration_miRNAs']) and not pd.isna(row['ML_max_pred_log2FC'])) else ('machine_learning' if pd.isna(row['manual_integration_miRNAs']) else 'manual_integration'), axis=1)
df['Gene name'] = df.index.map(gid2n)

writer = pd.ExcelWriter(snakemake.output['supp_table'], engine='xlsxwriter')
df.to_excel(writer, sheet_name='Main', startrow=1, index=True)
worksheet = writer.sheets['Main']

worksheet.write(0, 0, 'Supp. Table 8: miRNA target genes as predicted by filtering approach and/or by machine learning approach.')
writer.save()
