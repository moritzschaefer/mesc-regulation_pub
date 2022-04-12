import pandas as pd

from moritzsphd.external.heap_ago2 import bed_export_mre_mirnas

heap_ago2_sites = pd.read_csv(snakemake.input['ago2_heap_data'], index_col=0)

mirna_expr = pd.read_excel(snakemake.input['mirna_data'], skiprows=2, index_col=0)[['RIP_AGO2', 'RIP_AGO1']].mean(axis=1)
selector = (mirna_expr.reindex(heap_ago2_sites.mirna).fillna(0) > 10).values

bed_export_mre_mirnas(snakemake.output[0], heap_ago2_sites.loc[selector])
