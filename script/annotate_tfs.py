import pandas as pd

try:
    tf_annotation = snakemake.input['tf_annotation']
    data = snakemake.input['data']
    out = snakemake.output['data']
except NameError:
    tf_annotation = '../output/tf_annotation.csv'
    data = '../output/mirnas/interaction_ranking_protein_coding.csv'
    out = '../output/mirnas/interaction_ranking_protein_coding_istf.csv'

tfs = pd.read_csv(tf_annotation, index_col=0)
df = pd.read_csv(data, index_col=0)
df.insert(0, 'is_tf', df.index.map(lambda ensid: ensid in tfs.index))

df.to_csv(out)
