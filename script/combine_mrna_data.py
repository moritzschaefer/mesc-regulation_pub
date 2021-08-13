import itertools

import numpy as np
import pandas as pd


def cut_version(s: str):
    if type(s) is str:
        if '.' in s:
            return s[:s.find('.')]
        else:
            return s
    else:
        try:
            name = s.name
        except AttributeError:
            name = None
        return pd.Series([cut_version(v) for v in s], name=name)


mutants = snakemake.params['mutants']

tpm_expression = pd.read_csv(snakemake.input['tpm_expression'], sep='\t', index_col=0)

dfs = {}
for mutant, diffexp_f in zip(mutants, snakemake.input['diffexp']):
    dfs[mutant] = pd.read_csv(diffexp_f, sep='\t', index_col=0)
    dfs[mutant].index = cut_version(dfs[mutant].index)
    repls = [f'{mutant}_1', f'{mutant}_2']
    dfs[mutant]['tpm_expression'] = tpm_expression[repls].mean(axis=1)

repls = [f'WT_1', f'WT_2']
dfs['WT'] = pd.DataFrame({'log2FoldChange': 0.0,
                          'padj': 1.0,
                          'tpm_expression': tpm_expression[repls].mean(axis=1).reindex(dfs[mutants[0]].index)})

columns = list(itertools.product(snakemake.params['mutants'] + ['WT'], ['tpm_expression', 'log2FoldChange', 'padj']))
df = pd.DataFrame(
    index=pd.MultiIndex.from_arrays((dfs[mutants[0]].index, dfs[mutants[0]]['external_gene_name'].values)),
    columns=pd.MultiIndex.from_tuples(columns),
    data=np.array([dfs[mutant][field].values for mutant, field in columns]).T
)
df.index.name = 'Geneid'

df.to_csv(snakemake.output['csv'])

writer = pd.ExcelWriter(snakemake.output['supp_table'], engine='xlsxwriter')
df.to_excel(writer, sheet_name='Main', startrow=1, index=True)
worksheet = writer.sheets['Main']

worksheet.write(0, 0, 'Supp. Table 1: RNA-seq TPM values and differential expression for WT and/versus RNAi KO mutants. Referenced by (Supp.) Figures 1, 2 and 3')
writer.save()
