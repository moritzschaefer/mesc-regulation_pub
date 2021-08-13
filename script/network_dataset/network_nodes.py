import logging
import re

import gffutils
import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

sample = snakemake.wildcards['sample']
sample_1_2 = [f'{sample}_{i}' for i in [1, 2]]

db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)

# Load data
mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1]) \
              .reset_index(level=1, drop=True)
mirna_expression = pd.read_csv(snakemake.input['mirna_expression'], index_col=0, header=[0, 1]) \
    .xs('Expression', axis=1, level=1)
mirna_interactions = pd.read_csv(snakemake.input['mirna_edges'], index_col=[0, 1])
tf_interactions = pd.read_csv(snakemake.input['tf_edges'], index_col=[0, 1])

# Get gene names with WT expression
genes = tf_interactions.index.get_level_values(0).union(
    tf_interactions.index.get_level_values(1).union(
        mirna_interactions.index.get_level_values(1))).unique()

genes = mrna_data.loc[
    genes,
    (sample, 'tpm_expression')
].to_frame()
genes.columns = genes.columns.droplevel(0)

genes.rename(columns={'tpm_expression': 'wt_expression'}, inplace=True)

genes['wt_expression'] = MinMaxScaler().fit_transform(X=np.log2(genes['wt_expression'] + 1).values.reshape(-1, 1))
genes['is_gene'] = True
genes['is_tf'] = genes.index.map(lambda gid: gid in tf_interactions.index.get_level_values(0))
genes.index.name = 'node_id'

# Get miRNA names with WT expression
mirna_wt_expression = mirna_expression[sample]

mirnas = pd.DataFrame(index=mirna_interactions.index.get_level_values(0).unique())
mirnas['is_mirna'] = True
mirnas['wt_expression'] = mirna_wt_expression.loc[mirnas.index]
mirnas['wt_expression'] = MinMaxScaler().fit_transform(X=np.log2(mirnas['wt_expression'] + 1).values.reshape(-1, 1))
mirnas.index.name = 'node_id'

# Gene-region accessibilities

# In Bash it's as I would expect it to be
bedtools_output_names = [
    'chrom', 'source', 'feature', 'start', 'end', 'dummy1', 'strand', 'dummy2', 'attribs',
    'Ago12_1', 'Ago12_2', 'Ago1_1', 'Ago1_2', 'Ago2_1', 'Ago2_2', 'Dgcr8_1', 'Dgcr8_2',
    'Drosha_1', 'Drosha_2', 'Dicer_1', 'Dicer_2', 'WT_1', 'WT_2', 'WT2i_1', 'WT2i_2'
]

def _replace_mirna_ids(row):
    '''
    replace where mirbase with mirbase IDs
    as ensembl/gencode just has premature-miRNAs, I have to extrapolate to -3p and 5ps
    '''
    ret = [row['gene_id']]

    gn = row['gene_name']

    # the following miRNAs didn't get assigned properly (version issue..), so I do it manually:
    manual = {
        'Mir292': ['mmu-miR-292a-3p', 'mmu-miR-292a-5p'],
        'Mir92-1': ['mmu-miR-92a-3p'],
        'Mir301': ['mmu-miR-301a-3p', 'mmu-miR-301a-5p'],
        'Mir497': ['mmu-miR-497a-5p'],
        'Mir7-1': ['mmu-miR-7a-5p', 'mmu-miR-7a-1-3p'],
        'Mir124a-1': ['mmu-miR-124-3p', 'mmu-miR-124-5p'],
        'Mir142': ['mmu-miR-142a-5p', 'mmu-miR-142a-3p'],
        # ['mmu-miR-935'],  doesn't exist
        'Mir297-1': ['mmu-miR-297a-3p'],
        'Mir92-1': ['mmu-miR-92a-2-5p', 'mmu-miR-92a-1-5p'],
        'Mir3102': ['mmu-miR-3102-5p.2-5p', 'mmu-miR-3102-3p.2-3p'],
        'Mir466f-1': ['mmu-miR-466f'],
        # ['mmu-miR-9769-3p'], doesn't exist
        'Mir486': ['mmu-miR-486b-5p', 'mmu-miR-486a-5p'],
        'Mir18': ['mmu-miR-18a-5p', 'mmu-miR-18a-3p'],
        'Mir1191': ['mmu-miR-1191a'],
        'Mir465': ['mmu-miR-465a-3p'],
        'Mir450-1': ['mmu-miR-450a-5p']
    }

    if row['gene_source'] == 'mirbase' and not gn.startswith('Gm'):
        base = 'mmu-let-' if 'let' in gn else 'mmu-miR-'
        mirna = re.search('\d.*', gn).group()
        ret.extend([
            base + mirna + '-3p',
            base + mirna + '-5p',
        base + mirna])
        if mirna.endswith('-1'):
            ret.extend([
                base + mirna[:-2] + '-3p',
                base + mirna[:-2] + '-5p'
            ])
        try:
            ret.extend(manual[gn])
        except KeyError:
            pass

    return ret

def load_bedtools_output(fn):
    df = pd.read_csv(fn, header=None, names=bedtools_output_names, sep='\t', comment='#', index_col=False)

    attrs = df['attribs'].str.split('; ').apply(
        lambda row: {p[:p.find(' ')]: p[p.find(' ')+1:].strip('"')
                        for p in row})
    return df.join(pd.DataFrame(list(attrs)))


gene_acc_df['accessibility'] = gene_acc_df[sample_1_2].mean(axis=1)
# normalize by gene length
gene_acc_df['accessibility'] = gene_acc_df.apply(lambda v: v['accessibility'] / len(db[v['gene_id']]), axis=1)

gene_acc_df['node_id'] = gene_acc_df.apply(_replace_mirna_ids, axis=1)
gene_acc_df = gene_acc_df.explode('node_id')
gene_acc_df.drop(gene_acc_df.index[gene_acc_df.node_id.isna()], inplace=True)
gene_acc_df.set_index('node_id', inplace=True)

# TSS accessibilities
# taking the maximum among all transcripts of a gene makes most sense for TSS


tss_acc_df['accessibility'] = tss_acc_df[sample_1_2].mean(axis=1)
tss_acc_df = tss_acc_df.groupby(['gene_id', 'gene_name', 'gene_source'])['accessibility'].max().reset_index()

tss_acc_df['node_id'] = tss_acc_df.apply(_replace_mirna_ids, axis=1)
tss_acc_df = tss_acc_df.explode('node_id')
# tss_acc_df.drop(tss_acc_df.index[tss_acc_df.node_id.isna()], inplace=True)  # not necessary anymore
tss_acc_df.set_index('node_id', inplace=True)

# H3K27me3
h3k27me3_tss_df = pd.read_csv(snakemake.input['h3k27me3_tss_coverage'], index_col=0)


# combine output
nodes = pd.concat([genes, mirnas]).fillna(False)
nodes['gene_accessibility'] = MinMaxScaler().fit_transform(X=np.log2(gene_acc_df.reindex(nodes.index)['accessibility']+1).values.reshape(-1, 1))
nodes['tss_accessibility'] = MinMaxScaler().fit_transform(X=np.log2(tss_acc_df.reindex(nodes.index)['accessibility']+1).values.reshape(-1, 1))
# Many signals are negative (because of the log2 and IP-input being stronger than the IP itself)
# It therefore makes sense to apply 0-1 transformation
nodes['h3k27me3_tss'] = MinMaxScaler().fit_transform(X=h3k27me3_tss_df.reindex(nodes.index)[sample].values.reshape(-1, 1))

try:
    histone_df = pd.concat([pd.read_csv(h, index_col=0) for h in snakemake.input['external_histone_data']], axis=1)

    histone_df = pd.DataFrame(data=MinMaxScaler().fit_transform(X=np.log2(histone_df + 1)),
                              columns=snakemake.params['histone_marks'], index=histone_df.index)
    nodes = nodes.join(histone_df)

except ValueError:
    pass

logging.info(f'{nodes.isna().sum().sum()} nan values were found in network nodes and replaced with 0')
nodes.fillna(0, inplace=True)

nodes.to_csv(snakemake.output['nodes'])
