import pandas as pd

from scoring import score

try:
    df = pd.read_csv(snakemake.input['raw_interactions'], index_col=list(range(6)))
    # drop rows that don't have a genome position
    df = df.loc[df.index.get_level_values('start') >= 0]

    mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
    mrna_data = mrna_data[snakemake.params['mutants'] + ['WT']]
    log2fc = mrna_data.xs('log2FoldChange', axis=1, level=1)  # .loc[:, (slice(None), 'log2FoldChange')]
    padj = mrna_data.xs('padj', axis=1, level=1)
    tpm = mrna_data.xs('tpm_expression', axis=1, level=1)

    # TODO maybe I should use something simpler..
    up_genes = mrna_data.index[((log2fc > 0) &
                                (padj < snakemake.params['padj_threshold'])).sum(axis=1) >= snakemake.params['min_num_up_genes']]
    low_up_genes = mrna_data.index[(tpm > 1).any(axis=1) & (((log2fc > 0.1) & (log2fc < 0.5)).sum(axis=1) >= 4)]
    # low_up_genes = mrna_data.index[((log2fc > 0.1) & (log2fc < 0.5)).sum(axis=1) >= snakemake.params['min_num_up_genes']]

    down_genes = mrna_data.index[((log2fc < 0) &
                                  (padj < snakemake.params['padj_threshold'])).sum(axis=1) >= snakemake.params['min_num_up_genes']]

    with open(snakemake.output.up_genes, 'w') as f:
        f.write(','.join(up_genes.get_level_values(0)))
    with open(snakemake.output.low_up_genes, 'w') as f:
        f.write(','.join(low_up_genes.get_level_values(0)))
    with open(snakemake.output.down_genes, 'w') as f:
        f.write(','.join(down_genes.get_level_values(0)))

    for filtered, filename in zip([True, 'low', False], (snakemake.output['filtered'], snakemake.output['low_filtered'], snakemake.output['unfiltered'])):
        subdf = df.copy()

        if filtered:
            # filter for HEAP peaks (because we don't filter in combine_direct_mirna_interactions anymore..)
            subdf = subdf[subdf['AGO2 HEAP peak'] > 0]

            # filter for up-genes
            if filtered == 'low':
                subdf = subdf[subdf.index.get_level_values(0).isin(low_up_genes.get_level_values(0))]
            else:
                subdf = subdf[subdf.index.get_level_values(0).isin(up_genes.get_level_values(0))]

            # filter for miRNA expression
            subdf = subdf.loc[subdf['WT miRNA expression'] > snakemake.params['min_mirna_expression']]
            # filter for context++ score

            # this might not be ideal. Maybe there are collaborative interactions with individually low context++ scores
            subdf = subdf.loc[(subdf['weighted context++ score'] < snakemake.params['max_ts_score']) |
                              ((~subdf['is_3putr'].astype(bool)) & (subdf['MRE type'].isin(['7merm8', '8mer'])))]

        # Scorify (0.0-1.0) our features
        score_df, gene_order = score(subdf, padj, log2fc, snakemake.params['padj_threshold'])
        subdf = subdf.loc[score_df.index]
        subdf['Interaction score'] = score_df['interaction_score']
        subdf.reset_index().to_csv(filename, index=False)

        if filtered is True:
            # also store the by-gene-grouped ranking
            gene_order.index = pd.MultiIndex.from_arrays([gene_order.index, gene_order.index.map(lambda i: mrna_data.reset_index(level=1, drop=False).loc[i, 'level_1'][0])])
            gene_order.index.rename('Gene name', level=1, inplace=True)
            gene_order.name = 'score'

            gene_order.reset_index().to_csv(snakemake.output['grouped'], index=False)
        else:
            pass
            # TODO grouped unfiltered, but delete non-expressed miRNAs
except:
    import pdb; pdb.post_mortem()
    raise
