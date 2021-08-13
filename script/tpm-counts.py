'''
This code is an adaption from my previous rna-seq-star-deseq2 pipeline
'''
from pathlib import Path

import gffutils
import numpy as np
import pandas as pd

from moritzsphd.util.transcript_lengths import calculate_gene_lengths


def tpm(
        gene_counts,
        gene_lengths,
        mean_fragment_lengths=None,
):
    """Normalize the DGEList to transcripts per million.
    Adapted from Wagner, et al. 'Measurement of mRNA abundance using RNA-seq data:
    RPKM measure is inconsistent among samples.' doi:10.1007/s12064-012-0162-3
    Read counts :math:`X_i` (for each gene :math:`i` with gene length :math:`\widetilde{l_j}` )
    are normalized as follows:
    .. math::
        TPM_i = \\frac{X_i}{\\widetilde{l_i}}\cdot \\
        \\left(\\frac{1}{\sum_j \\frac{X_j}{\widetilde{l_j}}}\\right) \cdot 10^6
    Args:
        gene_lengths: 1D array of gene lengths for each gene.
        prior_count:
        mean_fragment_length: mean fragment length of reads
            (optional)
    """

    # compute effective length not allowing negative lengths
    index = gene_counts.index
    columns = gene_counts.columns
    if mean_fragment_lengths is None:
        effective_lengths = gene_lengths[:, np.newaxis]
    else:
        effective_lengths = (
            gene_lengths[:, np.newaxis] - mean_fragment_lengths[np.newaxis, :]
        ).clip(min=1)

    # how many counts per base
    base_counts = gene_counts / effective_lengths

    counts = 1e6 * base_counts.to_numpy() / np.sum(base_counts, axis=0)[np.newaxis, :]
    return pd.DataFrame(index=index, columns=columns, data=counts)


def extract_mean_fragment_lengths(filenames, samples):
    FIELD = 'MEAN_INSERT_SIZE'
    mean_fragment_lengths = pd.Series()
    for sample, filename in zip(samples,[Path(f) for f in filenames]):
        with filename.open() as f:
            index = None
            for l in f:
                if FIELD in l:
                    index = l.split('\t').index(FIELD)
                elif index is not None:
                    mean_fragment_lengths[sample] = float(l.split('\t')[index])
                    break

    return mean_fragment_lengths


counts = pd.read_csv(snakemake.input['expression'], sep='\t', index_col=0)

db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)
gene_lengths = calculate_gene_lengths(counts.index, db)
# mean_fragment_lengths = extract_mean_fragment_lengths(
#     snakemake.input.mean_frag_lengths,
#     snakemake.params['samples']
#     )

# we ignore mean_fragment_lengths because it leads to "unnaturally" high deviations.
tpm(counts, gene_lengths).to_csv(snakemake.output[0], sep="\t", index=True)
