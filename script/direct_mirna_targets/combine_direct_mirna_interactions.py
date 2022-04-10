try:
    import gffutils
    import pandas as pd
    import pyBigWig
    # import swifter
    from tqdm import tqdm
    from tqdm.contrib.concurrent import process_map
    tqdm.pandas()

    # load data
    mrna_data = pd.read_csv(snakemake.input['mrna_data'], index_col=[0, 1], header=[0, 1])
    mrna_data = mrna_data[snakemake.params['mutants'] + ['WT']]
    # mirna_expression = pd.read_csv(snakemake.input['mirna_expression'], index_col=0, header=[0, 1]) \
    #     .xs('Expression', axis=1, level=1)
    mirna_loading = pd.read_csv(snakemake.input['mirna_loading'], sep='\t', index_col=list(range(6)))
    mirna_loading.index = mirna_loading.index.get_level_values(0)
    mirna_loading = mirna_loading[['RIP_E14_AGO2', 'RIP_E14_AGO1']].mean(axis=1)
    ago2_heap_data = pd.read_csv(snakemake.input['ago2_heap_data'], index_col=0)  # of note, 6mers are not filtered here!

    # this leads to:  numpy/lib/arraysetops.py:580: FutureWarning: elementwise comparison failed; returning scalar instead,
    # but in the future will perform elementwise comparison
    ts_predictions = pd.read_csv(snakemake.input['targetscan_prediction'], index_col=0)
    # associate ts_predictions with heap coverage data

    # get up-genes
    log2fc = mrna_data.xs('log2FoldChange', axis=1, level=1)  # .loc[:, (slice(None), 'log2FoldChange')].
    padj = mrna_data.xs('padj', axis=1, level=1)
    expr = mrna_data.xs('tpm_expression', axis=1, level=1)
    log2fc.index = log2fc.index.droplevel(1)
    padj.index = padj.index.droplevel(1)
    expr.index = expr.index.droplevel(1)
    expressed_genes = set(expr.index[(expr > snakemake.params['min_mrna_expression']).any(axis=1)])

    # get potential mirna interactors
    # mirna_wt_expression = mirna_expression['WT']

    def _sum_scores(l, ascending=False):
        '''
        Aggregate multiple scores for the same gene.
        Avoid very high scores for targets with many targeting mirnas.
        Example:
        score = 1 * strongestinteraction + 0.5 second + 0.250* third
        '''
        return sum(l.sort_values(ascending=ascending) * [2**(-i) for i in range(len(l))])

    ts_predictions['is_3putr'] = True
    ts_predictions.rename(columns={'mre_start': 'start', 'mre_end': 'end'}, inplace=True)
    ts_predictions.set_index(['Geneid', 'miRNA', 'start', 'end'], inplace=True)

    # join HEAP peaks
    ago2_heap_data.rename(columns={'gene_id': 'Geneid', 'mirna': 'miRNA'}, inplace=True)
    ago2_heap_data = ago2_heap_data[ago2_heap_data.mre_type != '6mer']  # delete 6mers..
    ago2_heap_data.set_index(['Geneid', 'miRNA', 'start', 'end'], inplace=True)

    df = ts_predictions.join(ago2_heap_data, how='outer', lsuffix='_ts', rsuffix='_heap').fillna(0)

    # filter genes that are not in our data (like CDR1as from TargetScan)
    df = df.loc[df.index.get_level_values(0).isin(mrna_data.index.get_level_values(0))]

    # set chromosom for all entries (it was missing for the TargetScan ones)
    db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)
    gene_chromosomes = pd.Series({gene.id: gene.chrom for gene in db.features_of_type('gene')})
    df['chr'] = gene_chromosomes.reindex(df.index.get_level_values('Geneid')).values

    # An interaction is in the 3pUTR if it was detected by TS or by our own analysis (heap)
    df['is_3putr'] = df['is_3putr_ts'] | df['is_3putr_heap'].astype(bool)

    df['mre_type'] = df['mre_type_ts']
    df.loc[df['mre_type_ts'] == 0.0, 'mre_type'] = df.loc[df['mre_type_ts'] == 0.0, 'mre_type_heap']
    df.drop(columns=['is_3putr_ts', 'is_3putr_heap', 'mre_type_ts', 'mre_type_heap', 'Gene Symbol'], inplace=True)

    # # combine if a miRNA is having multiple seeds of the same type in the same gene_location for a given gene (shouldn't happen too often)
    # ago2_heap_data = ago2_heap_data.groupby(['Geneid', 'miRNA', 'mre_type', 'gene_location']) \
    #                                .agg({'score': _sum_scores, 'is_cds': 'max', 'is_5putr': 'max',
    #                                      'is_3putr': 'max', 'isoform_expression': 'max'})
    # # for simplicity, we neglect the few cases with multiple binding sites of the same type in a 3p UTR
    # ts_predictions = ts_predictions.groupby(['Geneid', 'miRNA', 'mre_type', 'gene_location']) \
    #                                .agg({'conserved': 'max', 'weighted context++ score': lambda v: _sum_scores(v, True)})

    # df['is_3putr'] = (df['is_3putr'].astype(bool) | (df.index.get_level_values('gene_location') == '3putr')).astype(int)
    df['gene_location'] = df.progress_apply(lambda row: '3putr' if row.is_3putr else
                                                        ('5putr' if row.is_5putr else 'cds'), axis=1)
    df.set_index('gene_location', append=True, inplace=True)

    # mutant_abbrv = {'Dicer': 'Di', 'Drosha': 'Dr', 'Dgcr8': 'Dg', 'Ago12': 'Ag'}
    df['upregulated mutants'] = ((padj.loc[df.index.get_level_values(0)] < snakemake.params['combined_padj_threshold']) &
                                (log2fc.loc[df.index.get_level_values(0)] > snakemake.params['log2fc_threshold'])).reset_index(drop=True).apply(  # convert to short-names
        lambda row: ','.join(row.index[row].map(lambda i: i[:2])),
        axis=1).values

    df['WT miRNA loading'] = mirna_loading.reindex(df.index.get_level_values(1)).values
    df['WT mRNA expression']  = expr.reindex(df.index.get_level_values(0))['WT'].fillna(0.0).values

    # Beautify data frame a bit before export
    df.rename(columns={'conserved': 'MRE conserved', 'score': 'AGO2 HEAP peak'}, inplace=True)
    df['Gene name'] = mrna_data.index.to_frame().reset_index(drop=True).set_index(0)[1].reindex(df.index.get_level_values(0)).values
    df.set_index('Gene name', append=True, inplace=True)

    df = df.reorder_levels(['Geneid', 'Gene name', 'gene_location', 'start', 'end', 'miRNA'])
    df.rename(columns={'mre_type': 'MRE type'}, inplace=True)

    df.reset_index().to_csv(snakemake.output['main'], index=False)
except:
    import pdb; pdb.post_mortem()
    raise
