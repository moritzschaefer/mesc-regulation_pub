import gffutils
import pandas as pd
from moritzsphd.data import average_sample
from tqdm import tqdm

tqdm.pandas()

try:
    df = pd.read_csv(snakemake.input['interactions'])
    # prepare isoform contributions from sequencing
    db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)

    rna_seq_ts_expr = pd.read_csv(snakemake.input['rna_seq_transcript_tpm'], sep='\t', index_col=0)
    quant_seq_ts_expr = pd.read_csv(snakemake.input['quant_seq_transcript_expr'], sep='\t', index_col=0)
    quant_seq_ts_expr *= (1e6 / quant_seq_ts_expr.sum())  # CPM-normalize the raw data

    rna_tdf = average_sample(rna_seq_ts_expr)[['WT']]
    rna_tdf['gene_id'] = rna_tdf.index.map(lambda tid: db[tid].attributes['gene_id'][0])
    rna_tdf_grouped = rna_tdf.groupby('gene_id')['WT'].apply(lambda g: g / g.sum()).to_frame()  # just normalize the sum
    rna_tdf_grouped['gene_id'] = rna_tdf_grouped.index.map(lambda tid: db[tid].attributes['gene_id'][0])

    quant_tdf = average_sample(quant_seq_ts_expr)[['WT']]
    quant_tdf['gene_id'] = quant_tdf.index.map(lambda tid: db[tid].attributes['gene_id'][0])
    quant_tdf_grouped = quant_tdf.groupby('gene_id')['WT'].apply(lambda g: g / g.sum()).to_frame()  # just normalize the sum
    quant_tdf_grouped['gene_id'] = quant_tdf_grouped.index.map(lambda tid: db[tid].attributes['gene_id'][0])

    def _isoform_expression(row):
        # find transcript contributions
        hits = db.region(seqid=str(row['chr']).replace('chr', ''), start=row['start'], end=row['end'])

        # only consider exons of given gene and take their transcript_ids
        tids = [hit.attributes['transcript_id'][0] for hit in hits if hit.featuretype == 'exon' and hit.attributes['gene_id'][0] == row['Geneid']]

        rna_iso = rna_tdf_grouped['WT'].reindex(tids).dropna().sum()
        quant_iso = quant_tdf_grouped['WT'].reindex(tids).dropna().sum()
        return rna_iso, quant_iso

    contribs = df.progress_apply(_isoform_expression, axis=1)
    df['rnaseq_isoform_expression'] = [c[0] for c in contribs]
    df['quantseq_isoform_expression'] = [c[1] for c in contribs]

    df.to_csv(snakemake.output['interactions'], index=False)
except:
    import pdb
    pdb.post_mortem()
    raise
