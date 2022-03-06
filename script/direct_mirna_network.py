import json

import networkx as nx
import numpy as np
import pandas as pd
from moritzsphd.data import mirbase_seqs
from scipy.stats.mstats import gmean

# TODO add information about PRC2, SwiSNF, ..., and other complexes
# TODO legend!

try:
    mutants = snakemake.params['mutants']
    interaction_fn = snakemake.input['interaction_fn']
    tf_fn = snakemake.input['tf_fn']
    mesc_mirna_fn = snakemake.input['mesc_mirna_fn']
    mrna_data_fn = snakemake.input['mrna_data_fn']
except NameError:
    mutants = ['Drosha', 'Dicer', 'Ago12']
    interaction_fn = '../output/mirnas/interaction_ranking_protein_coding.csv'
    tf_fn = '../output/tf_annotation.csv'
    mesc_mirna_fn = '../output/mesc_mirnas.csv'
    mrna_data_fn = '../output/mrna_data_protein_coding.csv'

df = pd.read_csv(interaction_fn)
tf_df = pd.read_csv(tf_fn, index_col=0)
mesc_mirnas = pd.read_csv(mesc_mirna_fn, index_col=0, squeeze=True)
mrna_data = pd.read_csv(mrna_data_fn, index_col=[0, 1], header=[0, 1])

# TODO ignore 2i for now
mrna_data = mrna_data[snakemake.params['mutants'] + ['WT']]

# invert p-value for downregulated genes (penalizing!)
for col in mrna_data.columns.get_level_values(0):
    downregulated = mrna_data.loc[:, (col, 'log2FoldChange')] < 0
    mrna_data.loc[downregulated, (col, 'padj')] = 1 / mrna_data.loc[downregulated, (col, 'padj')]
    # fix missing values:
    mrna_data.loc[:, (col, 'log2FoldChange')].fillna(0, inplace=True)
    mrna_data.loc[:, (col, 'padj')].fillna(1.0, inplace=True)

# mark ESC-specific miRNAs with dashed border (TODO also, mark ESC-specific TFs with dashed border)
df['mESC specific miRNA'] = mesc_mirnas.reindex(df['miRNA']).fillna(False).values

# filter nodes (i.e. combine miRNAs with same seed, just get max values when grouping and sum up the WT expression of the miRNAs)
def seed_for_mirna(mirna):
    mirna_seq = mirbase_seqs()[mirna]
    return str(mirna_seq[1:8])

# merge miRNA families
df['seed'] = df['miRNA'].apply(seed_for_mirna)
families = df.groupby('seed')['miRNA'].apply(lambda mirs: 'miR-' + '/'.join([mir.lstrip('mmu-').lstrip('-miR') for mir in np.unique(mirs)]))
df['miRNA family'] = families.loc[df['seed']].values
df['#(upregulated mutants)'] = df['upregulated mutants'].apply(lambda v: v.split().count(',') + 1)

df = df.groupby(['Geneid', 'Gene name', 'seed', 'miRNA family']).agg({
    'MRE conserved': 'max',
    'AGO2 HEAP peak': 'max',
    'MRE type': 'max',
    'gene_location': 'min',
    '#(upregulated mutants)': 'max',
    'Interaction score': 'max',
    'mESC specific miRNA': lambda group: 'yes' if np.mean(group) >= 0.5 else 'no',
    'weighted context++ score': 'min',
    'WT miRNA expression': 'sum'
})

# Upregulation Score (border size is upregulation-pvalue)

df['upregulation_score'] = -np.log10(mrna_data.loc[:, (slice(None), 'padj')].apply(gmean, axis=1))
df = df.reset_index()

# WT expression size is WT log2-expression (for both miRNAs and genes)
wt_tpm = mrna_data[('WT', 'tpm_expression')]
wt_tpm.index = wt_tpm.index.droplevel(1)
df['WT gene expression'] = np.log2(wt_tpm.loc[df['Geneid']] + 1).values
df['WT miRNA expression'] = np.log2(df['WT miRNA expression'] + 1)

# node type: TF or not
df = df
df['is_tf'] = 'no'
df.loc[df['Geneid'].isin(tf_df.index), 'is_tf'] = 'yes'
df.loc[df['Geneid'].isin(tf_df.index[tf_df.esc_specific]), 'is_tf'] = 'yes,esc_specific'

# fix data types for JSON export
df.loc[~df['MRE conserved'].isna(), 'MRE conserved'] = df.loc[~df['MRE conserved'].isna(), 'MRE conserved'].astype(float)
df['#(upregulated mutants)'] = df['#(upregulated mutants)'].astype(float)

# delete some low-quality connections
df = df[df['Interaction score'] > 1.2]

graph = nx.convert_matrix.from_pandas_edgelist(df.reset_index(), 'miRNA family', 'Gene name', ['Interaction score', 'AGO2 HEAP peak', 'MRE type'])  # 'MRE conserved', 'weighted context++ score',  <- has NaNs
for key in graph.nodes:
    try:
        values = df[df['Gene name'] == key].iloc[0][[
            '#(upregulated mutants)',
            'upregulation_score',
            'WT gene expression',
            'is_tf'
        ]]
        values['type'] = 'gene' if values['is_tf'] == 'no' else 'tf'
        values['mesc_specific'] = 'yes' if 'specific' in values['is_tf'] else 'no'
        values['expression'] = values['WT gene expression']
    except IndexError:
        values = df[df['miRNA family'] == key].iloc[0][[
            'WT miRNA expression',
            'mESC specific miRNA'
        ]]
        values['type'] = 'miRNA'
        values['expression'] = values['WT miRNA expression']
        values['mesc_specific'] = values['mESC specific miRNA']
    graph.nodes[key].update(**values.to_dict())

with open(snakemake.output['network'], 'w') as f:
    json.dump(nx.readwrite.json_graph.cytoscape_data(graph), f)
df.to_csv(snakemake.output['table'])
