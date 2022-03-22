'''
Number of identified targets when integrating different kinds of features. TODO: this is very similar to rank_direct_mirna_interactions

Only those genes that are expressed in WT or any of the RNAi mutants with at least 1 TPM are being considered.

TargetScan designed its context++ score only for 3' UTRs of mRNAs. To prevent the exclusion of potential CDS- or 5' UTR-interactions, we decided not to use the TargetScan scores for filtering interactions.
'''

import itertools
import pdb
# already filterted for minimum of 1 TPM in any mutant
from collections import defaultdict
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import upsetplot

from scoring import score

try:
    df = pd.read_csv(snakemake.input['raw_interactions'], index_col=list(range(6)))

    mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
    log2fc = mrna_data.xs('log2FoldChange', axis=1, level=1)  # .loc[:, (slice(None), 'log2FoldChange')]
    padj = mrna_data.xs('padj', axis=1, level=1)

    with open(snakemake.input.up_genes, 'r') as f:
        up_genes = dict((g, True) for g in f.read().split(','))

    filters = {
        'AGO2 binding': lambda row: row['AGO2 HEAP peak'] > 0,
        'WT miRNA expression': lambda row: row['WT miRNA expression'] > snakemake.params['min_mirna_expression'],
        'Mutant upregulation': lambda row: up_genes.get(row.name[0], False),
        'TargetScan score': lambda row: (row['weighted context++ score'] < snakemake.params['max_ts_score']) or
                              ((not row['is_3putr']) and (row['MRE type'] in ['7merm8', '8mer']))
    }
    selectors = {f: df.apply(filters[f], axis=1) for f in filters}

    index = pd.DataFrame(
        columns=selectors.keys(),
        data=itertools.product(*([(0, 1)] * len(selectors)))
    )

    filtered_dfs = {}
    def _apply_filters(binary):
        combined_selector = pd.Series(data=True, index=df.index)
        for f in binary.index[binary.astype(bool)]:  # only consider the filters that are true
            combined_selector &= selectors[f]
        subdf = df.loc[combined_selector].copy()
        if sum(binary) >= 3:
            # Scorify (0.0-1.0) our features
            score_df, _ = score(subdf, padj, log2fc, snakemake.params['combined_padj_threshold'])
            subdf = subdf.loc[score_df.index]
            subdf['Interaction score'] = score_df['interaction_score']

            # rank the df
            filtered_dfs['_'.join(binary.index[binary == 1])] = subdf

        return len(subdf.index.get_level_values(0).drop_duplicates())

    data = index.apply(_apply_filters, axis=1)
    num_genes = pd.Series(index=pd.MultiIndex.from_frame(index),
                          data=data.values)
    num_genes.to_csv(snakemake.output['counts'])
    # generate excel sheet, put number of genes, number of interactions and wich filtering was used in the headline

    writer = pd.ExcelWriter(snakemake.output['interactions'], engine='xlsxwriter')
    for i, (filter_name, filtered_df) in enumerate(filtered_dfs.items()):
        filtered_df.reset_index().set_index(['Geneid', 'Gene name', 'miRNA'])[[
            'weighted context++ score', 'gene_location', 'AGO2 HEAP peak',
             'upregulated mutants', 'WT miRNA expression', 'WT mRNA expression', 'Interaction score'
        ]].to_excel(writer, sheet_name=f'Filtering {i+1}', startrow=1, index=True)

        filters = filter_name.split('_')
        writer.sheets[f'Filtering {i+1}'].write(0, 0, f'Supp. Table 3: Filtered interactions from integrative analysis. This sheet\'s filters: ({", ".join(filters)}), leading to {len(filtered_df)} interactions on {len(filtered_df.index.get_level_values(0).drop_duplicates())} target genes. Referenced by Figure 2')
    writer.save()

except:
    pdb.post_mortem()
    raise
