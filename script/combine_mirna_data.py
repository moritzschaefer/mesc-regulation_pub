import itertools

import numpy as np
import pandas as pd
from moritzsphd.data import mirna_diffexp, special_normalized_mirnas

mutants = snakemake.params['mutants']

expr = special_normalized_mirnas()

dfs = {}
for mutant in mutants:
    dfs[mutant] = mirna_diffexp(mutant)
    dfs[mutant]['Expression'] = expr[[f'{mutant}_1', f'{mutant}_2']].mean(axis=1)

dfs['WT'] = pd.DataFrame({'log2FoldChange': 0.0,
                          'padj': 1.0,
                          'Expression': expr[[f'WT_1', f'WT_2']]
                            .mean(axis=1)
                            .reindex(dfs[mutants[0]].index)})

columns = list(itertools.product(snakemake.params['mutants'] + ['WT'], ['Expression', 'log2FoldChange', 'padj']))
df = pd.DataFrame(
    index=dfs[mutants[0]].index,
    columns=pd.MultiIndex.from_tuples(columns),
    data=np.array([dfs[mutant][field].values for mutant, field in columns]).T
)
df.index.name = 'miRNA'
df.to_csv(snakemake.output['csv'])

writer = pd.ExcelWriter(snakemake.output['supp_table'], engine='xlsxwriter')
df.to_excel(writer, sheet_name='Main', startrow=1, index=True)
worksheet = writer.sheets['Main']

worksheet.write(0, 0, 'Supp. Table2: sRNA-seq (microRNAs) CPM values and differential expression for WT and/versus RNAi KO mutants. Referenced by (Supp.) Figures 1, 2 and 3')
writer.save()
