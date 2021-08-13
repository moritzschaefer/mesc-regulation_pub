'''
score miRNA interactions on a gene-wise basis
'''
import numpy as np
import pandas as pd
from scipy.stats.mstats import gmean
from sklearn.preprocessing import MinMaxScaler


def score(df, padj, log2fc, padj_threshold=0.1):
    score_df = pd.DataFrame(index=df.index)
    # TODO maybe take score from  script/mir290_ml/dataset.py
    score_df['mutant_up_score'] = (
        (
            (padj.loc[df.index.droplevel([2, 3, 4])] < padj_threshold) &
            (log2fc.loc[df.index.droplevel([2, 3, 4])] > 0)
        ).sum(axis=1)).values / 4

    score_df['wt_mirna_expression'] = MinMaxScaler().fit_transform(X=np.log2(df['WT miRNA expression'] + 1).values.reshape(-1, 1))
    # MRE score is a bit more complicated because of the noncoding genes (they dont have a context++-score)
    cs = df['context++ score']
    score_df['mre_score'] = MinMaxScaler().fit_transform(-cs.values.reshape(-1, 1))
    score_df['heap_ago2'] = MinMaxScaler().fit_transform(np.log2(df['AGO2 HEAP peak'] + 1).values.reshape(-1, 1))

    # Combine scores, order the targets. Multiple binding miRNAs are geomean-summed (see below) to reduce the power of multi-mirna targets
    score_df['mirna_combined'] = score_df[['wt_mirna_expression', 'mre_score', 'heap_ago2']].sum(axis=1)
    score_df['interaction_score'] = score_df[['wt_mirna_expression', 'mre_score', 'heap_ago2', 'mutant_up_score']].sum(axis=1)
    # filter: if no/low TS prediction (e.g. no 3'UTR MRE), then require a min-score
    # if filtered:
    #     score_df = score_df.loc[(score_df['interaction_score'] > snakemake.params['cds_min_score']) | (cs < snakemake.params['max_ts_score'])]

    # now sort such that top target genes come first
    score_df.reset_index(inplace=True)

    gene_order = score_df.groupby('Geneid').apply(lambda group: group['mutant_up_score'].mean()  + gmean([group['mirna_combined'].sum(), group['mirna_combined'].max()])).sort_values(ascending=False)
    score_df['Geneid'] = pd.Categorical(score_df['Geneid'], categories=gene_order.index[::-1], ordered=True)
    score_df = score_df.sort_values(['Geneid', 'interaction_score'], ascending=False).set_index(df.index.names)
    score_df.drop(score_df.index[score_df.index.duplicated()], inplace=True) # 48 unimportant duplicates

    return score_df, gene_order
