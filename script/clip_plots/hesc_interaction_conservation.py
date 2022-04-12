import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig
import pyensembl
import seaborn as sns
from moritzsphd.data import mirbase_seqs
from scipy.stats import pearsonr

bw = pyBigWig.open(snakemake.input['hesc_clip_bw'])

# Take our mirna interactions
df = pd.read_csv(snakemake.input['interaction_ranking'])

def _get_seed_match(row):
    if row['MRE type'] == '8mer':
        seed_match = 'T' + mirbase_seqs()[row['miRNA']][1:8].back_transcribe()
    elif row['MRE type'] == '7merm8':
        seed_match = mirbase_seqs()[row['miRNA']][1:8].back_transcribe()
    elif row['MRE type'] == '7merA1':
        seed_match = 'T' + mirbase_seqs()[row['miRNA']][1:7].back_transcribe()
    else:
        raise ValueError(row['MRE type'])

    seed_match = seed_match.reverse_complement()
    return str(seed_match)


def _get_seed(row):
    if row['MRE type'] == '8mer':
        seed_match = mirbase_seqs()[row['miRNA']][1:8].back_transcribe()
    elif row['MRE type'] == '7merm8':
        seed_match = mirbase_seqs()[row['miRNA']][1:8].back_transcribe()
    elif row['MRE type'] == '7merA1':
        seed_match = mirbase_seqs()[row['miRNA']][1:7].back_transcribe()
    else:
        raise ValueError(row['MRE type'])

    return str(seed_match)


df['seed_match'] = df.apply(_get_seed_match, axis=1)
df['seed'] = df.apply(_get_seed, axis=1)

human_ensembl = pyensembl.EnsemblRelease(93, pyensembl.species.human)

# def
# human_gene = human_ensembl.gene_by_id(ortholog_df.iloc[10000, 3])
# tr = human_gene.transcripts[0]
# human_gene.exons
# tr.coding_sequence_position_ranges
# tr.three_prime_utr_sequence
# three_prime_utr_start = tr.last_stop_codon_spliced_offset + 1
# exon = tr.exons[0]


def peak_size(gene_id, seed_match):
    peak_sizes = [0]
    try:
        human_gene = human_ensembl.gene_by_id(gene_id)
    except:
        return None
    for transcript in human_gene.transcripts:
        try:
            three_prime_utr_start = transcript.last_stop_codon_spliced_offset + 1
        except ValueError:
            continue
        for match in re.finditer(seed_match, transcript.three_prime_utr_sequence):
            # convert to genome coordinates
            covered = 0
            for exon in transcript.exons:
                relative_start = match.start() + three_prime_utr_start - covered
                if relative_start > exon.length:
                    covered += exon.length
                else:
                    if transcript.strand == '+':
                        human_genome_start = exon.start + relative_start
                    else:
                        human_genome_start = exon.end - relative_start - len(seed_match)
                    break
            else:
                raise ValueError('cant happen...')
            # if gene_id == 'ENSG00000090447':
            #     print(seed_match, human_genome_start)

            # get peak size
            try:
                params = f'chr{human_gene.contig}', human_genome_start, human_genome_start + len(seed_match)
                peak_sizes.append(bw.stats(*params)[0])
            except RuntimeError:
                raise

    return max(peak_sizes)


ortholog_df = pd.read_csv(snakemake.input['orthologs'], sep='\t')
merged = df.reset_index().merge(ortholog_df[['Human gene stable ID', 'Gene stable ID']].rename(columns={'Gene stable ID': 'Geneid'}), on='Geneid', how='inner')
peak_sizes = merged.apply(lambda row: peak_size(row['Human gene stable ID'], row['seed_match']), axis=1)
df['human_peak_size'] = peak_sizes.groupby(merged['index']).max()

hesc_mirna_expr = pd.read_excel(snakemake.input['hesc_mirna_expr'], index_col=0, skiprows=1).dropna(how='all')
hesc_mirna_expr['seed'] = hesc_mirna_expr['Mature Sequence'].apply(lambda v: v[1:7].replace('U', 'T'))
hesc_mirna_expr['seed_m8'] = hesc_mirna_expr['Mature Sequence'].apply(lambda v: v[1:8].replace('U', 'T'))

seed_cpms = pd.concat([hesc_mirna_expr.groupby('seed')['RPM ESC'].sum(),
                       hesc_mirna_expr.groupby('seed_m8')['RPM ESC'].sum()])

df['human_mirna_family_cpm'] = seed_cpms.reindex(df['seed'].values).fillna(0).values
df['log_human_peak_size'] = np.log2(df['human_peak_size'] + 1)
non_expressed_mirna_peaks = df.loc[df['human_mirna_family_cpm'] == 0, 'log_human_peak_size']
expressed_mirna_peaks = df.loc[df['human_mirna_family_cpm'] > 10, 'log_human_peak_size']
# enrichment-control
fig, ax = plt.subplots(figsize=(1, 4))

pd.DataFrame({'non-expressed human miRNAs': 100 * (non_expressed_mirna_peaks > 0).value_counts() / len(non_expressed_mirna_peaks),
                'expressed human miRNAs': 100 *(expressed_mirna_peaks > 0).value_counts() / len(expressed_mirna_peaks)}).loc[True].plot.bar(color='black')
ax.set_ylabel('% of human-conserved interactions with CLIP-seq reads')
fig.savefig(snakemake.output['mirna_expr_clip_conservation'])

# correlation:
fig, ax = plt.subplots(figsize=(3, 4))

pos_peaks = expressed_mirna_peaks = df.loc[(df['human_mirna_family_cpm'] > 10) & (df['human_peak_size'] > 0), ['human_mirna_family_cpm', 'human_peak_size']]
sns.regplot(data=np.log10(pos_peaks), x='human_mirna_family_cpm', y='human_peak_size', ax=ax, scatter_kws={'s': 5}, color='black')
sns.despine()
# ax.set_yticks([10, 100, 1000, 10000])
ax.set_xticklabels([f'{int(10**x)}' for x in ax.get_xticks()])
ax.set_yticks([0, 1, 2, 3])
ax.set_yticklabels([f'{int(10**x)}' for x in ax.get_yticks()])
print(f"Pearson correlation coefficient: {pearsonr((pos_peaks)['human_mirna_family_cpm'], (pos_peaks)['human_peak_size'])[0]:.3}")  # similar value for np.log2..
fig.savefig(snakemake.output['mirna_expr_clip_correlation'])
