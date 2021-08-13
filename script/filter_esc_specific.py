import pandas as pd

# input

mescspecific='output/mirnas/interaction_ranking_{subset}_mescspecific.csv',
nonspecific='output/mirnas/interaction_ranking_{subset}_nonspecific.csv',
mesc_mirnas = pd.read_csv(snakemake.input['mesc_mirnas'], index_col=0, squeeze=True)

df = pd.read_csv(snakemake.input['mirna_interaction_ranking'])

mescspecific_mask = df['miRNA'].isin(mesc_mirnas.index[mesc_mirnas.values].values)

df[mescspecific_mask].to_csv(snakemake.output['mescspecific'], index=False)
df[~mescspecific_mask].to_csv(snakemake.output['nonspecific'], index=False)
