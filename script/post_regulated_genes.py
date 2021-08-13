import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn3

from moritzsphd.data import ensembl_release


def _gene_from_transcript(tid):
    try:
        return ensembl_release.transcript_by_id(tid).gene_id
    except ValueError:
        return None

df = pd.read_csv(snakemake.input.degradation, index_col=0)
tf_df = pd.read_csv(snakemake.input.tf_annotation, index_col=[0, 1])
mrna_data = pd.read_csv(snakemake.input.mrna_data, index_col=[0, 1], header=[0, 1])
expr = mrna_data.xs('tpm_expression', axis=1, level=1)
expr.index = expr.index.droplevel(1)
expressed_genes = set(expr.index[(expr > 1).any(axis=1)])  # same filter as used in the RNAi_KO-based predictons..
print(len(expressed_genes))

mirna_regulated_transcripts = df[(df['padj'] < 0.5) & (df['log2fc'] > 0.5)]  # log2fc is inverted, > 0.5 means decreased in Dicer (maybe double check in raw data).
mirna_regulated_transcripts['Geneid'] = mirna_regulated_transcripts.index.map(_gene_from_transcript)
# only retain transcripts with associated gene
mirna_regulated_transcripts = mirna_regulated_transcripts.loc[~mirna_regulated_transcripts['Geneid'].isna()]  

pipeline_mirna_targets = pd.read_csv(snakemake.input.rna_based_predictions, index_col=0).query('score > 1')

venn3((set(pipeline_mirna_targets.index),  # already filtered for expressed genes..
       set(mirna_regulated_transcripts['Geneid']) & expressed_genes,
       set(tf_df.index.get_level_values(0)) & expressed_genes),
      ('RNAi-mutant derived miRNA targets (my pipeline)',
       '4sU derived miRNA targets (using log2fc filtering)',
       'Transcription Factors'))

plt.savefig(snakemake.output.overlap_plot)
mirna_regulated_transcripts.to_csv(snakemake.output.targets)
