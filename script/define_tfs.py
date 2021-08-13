import tarfile

import pandas as pd

mrna_df = pd.read_csv(snakemake.input['mrna_data'], header=[0, 1], index_col=[0, 1])
expressions = mrna_df.xs('tpm_expression', axis=1, level=1)

with tarfile.open(snakemake.input['tf_data'][0], "r:*") as tar:
    df = pd.read_csv(tar.extractfile('TF.expression.FPKM'), index_col=0, sep="\t")

tf_index = df.index.append(pd.Index(snakemake.params.add_tfs)).drop(pd.Index(snakemake.params.remove_tfs))

tf_max_expressions = expressions.max(axis=1).loc[(slice(None), tf_index)]
tfs = tf_max_expressions.index[tf_max_expressions > 2]  # at least 2 TPM in any of the mutants (~800 TFs)

# their ESC correlates very badly with our ESC mRNA seq counts (maybe because of TPM normalization?)
df_tfs = df.reindex(tfs.get_level_values(1))
esc_specific = (df_tfs['ESC'] == df_tfs.max(axis=1)) & \
   (df_tfs['ESC'] > df_tfs.mean(axis=1) * 2)

pd.DataFrame(index=tfs, data={
    'wt_expression': expressions.loc[tfs, 'WT'],
    'esc_specific': esc_specific.values,
}).to_csv(snakemake.output[0])
