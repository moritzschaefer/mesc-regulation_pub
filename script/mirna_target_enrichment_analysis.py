import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from moritzsphd.analysis.enrich import enrichr
from moritzsphd.data import gid2n

gene_sets = dict(zip(snakemake.params['gene_set_labels'], snakemake.input['gene_sets']))
enrichments = dict(zip(snakemake.params['gene_set_labels'], snakemake.input['enrichments']))
libraries = dict(zip(snakemake.params['gsea_libraries'], snakemake.input['libraries']))

human_mouse_trans = pd.read_csv(snakemake.input.human_mouse_trans, sep='\t', index_col=1).iloc[:, 0].to_dict()

gene_set_library = {}

for key, enrichment_file in enrichments.items():
    df = pd.read_csv(enrichment_file, index_col=0)
    for i, l in df.iterrows():
        with open(libraries[l.Gene_set]) as f:
            for line in f:
                fields = line.rstrip('\n').split('\t')
                if fields[0] == l.Term:
                    gene_set_library[f'{l.Gene_set}__{l.Term}'] = [human_mouse_trans[g] for g in fields[2:-1] if g in human_mouse_trans]

dfs = []
for key, fn in gene_sets.items():
    with open(fn) as f:
        genes = [gid2n(gid) for gid in f.read().rstrip('\n').split(',')]
    df, ax = enrichr(genes, gene_sets=gene_set_library, filter_padj=None, adjust_texts=True, plot=False, background='mmusculus_gene_ensembl')
    df['set'] = pd.Categorical([key] * len(df), categories=snakemake.params['gene_set_labels'], ordered=True)
    dfs.append(df)
df = pd.concat(dfs)
TOP_TERMS = 10

selected_terms = []
for s in df['set'].drop_duplicates():
    selected_terms.extend(df.loc[(df['set'] == s)].sort_values('-log10(adjusted p-value)', ascending=False).iloc[:TOP_TERMS]['Term'])

# put developmental biology first (it is among the top 8 hits of miRNA target genes)
selected_terms.remove('BioPlanet_2019__Developmental biology')
selected_terms.insert(0, 'BioPlanet_2019__Developmental biology')

df = df[df['Term'].isin(selected_terms)]
df['Term'] = pd.Categorical(df['Term'], categories=pd.unique(selected_terms), ordered=True)
df = df.sort_values('Term')
df['library'] = df['Term'].apply(lambda s: s.split('__')[0])
df['Term'] = df['Term'].apply(lambda s: re.sub(' \(.*\)', '', s.split('__')[1]).replace('positive', 'pos.').replace('regulation', 'reg.'))

fig, ax = plt.subplots(figsize=(5.1, 1))
df['resized_adjp'] = np.log2(df['-log10(adjusted p-value)'] + 1)  # NOTE: This requires manual fixing of the legend labels!!

ax = sns.scatterplot(data=df, x='Term', y='set', hue='library', size='resized_adjp', ax=ax, legend='brief', sizes=(0, 200), palette='muted')
ax.grid(False)
ax.set_ylim([-0.5, 2.5])
plt.tight_layout()
fig.savefig(snakemake.output.plot)
df.to_excel(snakemake.output.data)
