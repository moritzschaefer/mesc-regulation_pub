import itertools

import numpy as np
import pandas as pd

expr = pd.read_csv(snakemake.input['expr'], sep='\t', index_col=0)
samples = snakemake.params['samples']
neg = snakemake.params['neg_control']

matrix_cols = []
cpm = expr * (1e6 / expr.sum())  # CPM normalization (is the same as TPM for quant-seq)
dfs = {}
for sample, diffexp_f in zip(samples, snakemake.input['diffexps']):
    if 'miR' in sample:
        shortname = sample.replace('-', '')[:6] if 'miR' in sample else sample
        repls = [f'{shortname}KO_1', f'{shortname}KO_2']
    else:
        repls = [f'{sample}_1', f'{sample}_2']

    matrix_cols.extend(repls)

    dfs[sample] = pd.read_csv(diffexp_f, sep='\t', index_col=0)
    dfs[sample]['tpm_expression'] = cpm[repls].mean(axis=1)

dfs[neg] = pd.DataFrame({'log2FoldChange': 0.0,
                         'baseMean': -1.0,
                         'padj': 1.0,
                         'tpm_expression': cpm[[f'{neg}_1', f'{neg}_2']].mean(axis=1).reindex(dfs[samples[0]].index)})
matrix_cols.extend([f'{neg}_1', f'{neg}_2'])

# save expression matrix
expr[matrix_cols].to_csv(snakemake.output['read_count_matrix'], sep='\t')

columns = list(itertools.product(samples + [neg], ['tpm_expression', 'baseMean', 'log2FoldChange', 'padj']))
df = pd.DataFrame(
    index=pd.MultiIndex.from_arrays((dfs[samples[0]].index, dfs[samples[0]]['external_gene_name'].values)),
    columns=pd.MultiIndex.from_tuples(columns),
    data=np.array([dfs[sample][field].values for sample, field in columns]).T
)
df.index.name = 'Geneid'

writer = pd.ExcelWriter(snakemake.output['supp_table'], engine='xlsxwriter')
df.to_excel(writer, sheet_name='Main', startrow=1)
worksheet = writer.sheets['Main']

worksheet.write(0, 0, snakemake.params['title'])
writer.save()
