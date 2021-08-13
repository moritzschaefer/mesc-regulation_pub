import pandas as pd
from sklearn.preprocessing import MinMaxScaler

from moritzsphd.data import gn2id

# TODO needs refactoring since 'repressors' doesnt exist anymore. use what's defined in config.yaml (or better in selected_parameters)

# TODO some TFs have replicates. Their results should be combined somehow instead of just adding them twice!!
dfs = []
for tf, target_fn in zip(snakemake.params['tfs'].keys(), snakemake.input['tf_targets']):
    df = pd.read_csv(target_fn, index_col=0)
    df['tf'] = gn2id(tf.split('_')[0])
    # qvalue is stored in -log10-units.
    df['interaction_potential'] = MinMaxScaler().fit_transform(df['qvalue'].values.reshape(-1, 1))

    # interaction_potential has to have the sign of the regulation type
    if tf in snakemake.params['repressors']:
        df['interaction_potential'] *= -1
    dfs.append(df)

df = pd.concat(dfs)
df.rename(columns={'gene_id': 'target'}, inplace=True)

# TODO store distance/tss_distance

# the name 'interaction_potential' is concordant with what is used in the mirna network
df['interaction_type'] = 'tf'
df = df.set_index(['tf', 'target'])[['interaction_potential', 'interaction_type']]
df.to_csv(snakemake.output['edges'])
