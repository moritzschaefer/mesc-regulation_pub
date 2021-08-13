import itertools

import gffutils
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import hsv_to_rgb, rgb_to_hsv
from pybedtools import BedTool
from tqdm.contrib.concurrent import process_map

from moritzsphd.data import grc_seqs, mirbase_seqs
from moritzsphd.external.heap_ago2 import (bed_export_mre_mirnas,
                                           iterate_mirna_seqs, mre_target,
                                           peak_mirnas)

grc_seqs()  # preload chromosomal sequences
heap_df = pd.read_csv(snakemake.input['heap_peaks'])
heap_df['seqnames'] = heap_df['seqnames'].str.replace('chr', '')

mirna_df = pd.read_csv(snakemake.input['mirna_data'], header=[0, 1], index_col=0)
mirna_expr = mirna_df[('WT', 'Expression')]
mirnas = mirna_expr.index[mirna_expr > snakemake.params['min_mirna_expression']]

db = gffutils.FeatureDB(snakemake.input.annotation, keep_order=True)

# add missing fields in search regions
regions = snakemake.params['regions'].values()
for region in regions:
    region.update({'seqnames': region['chr'], 'score': '-1'})

searches = [(mirna, mirbase_seqs()[mirna], region) for mirna, region in
            itertools.product(mirnas, regions)]

binding_sites = process_map(iterate_mirna_seqs, searches)
binding_sites = [bs for sublist in binding_sites for bs in sublist if bs['mre_type'] != '6mer']

results = process_map(mre_target, binding_sites, chunksize=100)
results = [line for subdf in results for index, line in subdf.iterrows()]
df = pd.concat(results, axis=1, keys=[s.name for s in results], names=['gene_id']).T.reset_index()

# add missing zic2 peak (see [[file:~/wiki/gtd/reviews.org :ID:       4e989730-755c-4d9b-820b-f38651fe2f1f]])
zic2_mres = peak_mirnas({'seqnames': '14', 'start': 122479850, 'end': 122480330, 'strand': '.', 'score': -1},
                        min_mirna_expression=snakemake.params['min_mirna_expression'])
zic2_df = pd.DataFrame(zic2_mres)
zic2_df['gene_id'] = 'ENSMUSG00000061524'
zic2_df['is_cds'] = False
zic2_df['is_5putr'] = False
zic2_df['is_3putr'] = True
zic2_df = zic2_df.loc[zic2_df.mre_type != '6mer']


df = pd.concat((df, zic2_df)).reset_index(drop=True)
df['strand'] = df['mre_type'].apply({'7merA1': '+', '8mer': '.', '7merm8': '-'}.get)
df['score'] = df['mirna'].apply(mirna_expr.get)
df = df.groupby(['chr', 'start', 'end', 'mre_type', 'strand', 'offset']).agg({
    'score': np.sum,
    'mirna': lambda mres: 'miR-' + '/'.join([mre.lstrip('mmu-').lstrip('-miR') for mre in mres]),
}).reset_index()

df['score'] = np.log2(df['score'])

cmap = sns.color_palette("flare", as_cmap=True)

min_color = np.log2(snakemake.params['min_mirna_expression'])
max_color = np.log2(mirna_expr.max())
# BedTool.from
mre_bt = BedTool.from_dataframe(df.reset_index()[['chr', 'start', 'end', 'index']])
heap_bt = BedTool.from_dataframe(heap_df[['seqnames', 'start', 'end', 'name']])
in_heap_ids = mre_bt.intersect(heap_bt).to_dataframe()['name'].values
df['in_heap'] = False
df.loc[in_heap_ids, 'in_heap'] = True


def _color(row):
    col = [((row['score'] - min_color) / (max_color - min_color))] * 3
    # col = cmap((row['score'] - min_color) / (max_color - min_color))
    # col = col[:3]
    # if not row['in_heap']:
    #     hsv = rgb_to_hsv(col)
    #     hsv[2] = min(0.35 + hsv[2], 1.0)  # make it lighter to indicate it's not within a HEAP peak
    #     hsv[2] = min(0.35 + hsv[2], 1.0)  # make it lighter to indicate it's not within a HEAP peak
    #     col = hsv_to_rgb(hsv)

    return ','.join([str(int(v*255)) for v in col])


df['color'] = df.apply(_color, axis=1)

bed_export_mre_mirnas(snakemake.output[0], df)
