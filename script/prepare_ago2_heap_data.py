from moritzsphd.external.heap_ago2 import HeapAgo2

# get the binding sites

heap_ago2_sites = HeapAgo2(0, num_threads=snakemake.threads).mre_mirnas

# bools are annoying
for c in ['is_3putr', 'is_5putr', 'is_cds']:
    heap_ago2_sites[c] = heap_ago2_sites[c].fillna(0).astype(int)

# Zic2-annotation-fix
heap_ago2_sites.loc[heap_ago2_sites.gene_id == 'ENSMUSG00000061524', 'is_3putr'] = 1
heap_ago2_sites.to_csv(snakemake.output[0])
