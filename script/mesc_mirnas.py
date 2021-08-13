import re

import pandas as pd

sample_info = pd.read_csv(snakemake.input['sample_info'][0], sep='\t', index_col=0)

esc_samples = sample_info['description'].str.contains('Embryonic Stem cell').values
ips_samples = sample_info['description'].str.contains('iPS').values

mirna_cpms = pd.read_csv(snakemake.input['mirna_cpms'][0], sep='\t', index_col=0)

esc_mean = mirna_cpms.loc[:, esc_samples].mean(axis=1)
other_mean = mirna_cpms.loc[:, ~(esc_samples | ips_samples)].mean(axis=1)  # ignore ips_samples
other_max = mirna_cpms.loc[:, ~(esc_samples | ips_samples)].max(axis=1)  # ignore ips_samples

esc_mirna_selector = (esc_mean > (5 * other_mean))  # 5x works okish according to log output (see bottom of script)
esc_mirnas = mirna_cpms.index[esc_mirna_selector]
nonesc_mirnas = mirna_cpms.index[~esc_mirna_selector]

esc_mirna_selector.name = 'ESC specific'
esc_mirna_selector.index.name = 'miRNA'
esc_mirna_selector.to_csv(snakemake.output[0])

# TODO these are control mirnas from publication PMC4118578, I'll use them to control my analysis
control_mesc_mirnas = [
    'mir-290',
    'mir-291',
    'mir-292',
    'mir-293',
    'mir-294',
    'mir-295',
    'mir-296',
    'mir-302',
    'mir-17',
    'mir-92',
    'mir-15b',
    'mir-16',
    'mir-17',
]

def _mirna_in_control(mirna):
    for cmirna in control_mesc_mirnas:
        if re.match(f'({cmirna}$)|({cmirna}[^0-9])', mirna.lower()[4:]):
            return True
    return False

matched = esc_mirna_selector.index.map(_mirna_in_control).to_series(index=esc_mirna_selector.index)
with open(snakemake.log[0], 'w') as f:
    f.write(f'{matched.sum()} miRNAs identified as control. {(matched & esc_mirna_selector).sum()} from those were identified as mESC-specific.\n')
    f.write(f'Of the other {(matched == False).sum()} miRNAs, {((matched == False) & esc_mirna_selector).sum()} were positively identified as mESC-specific.\n')
    f.write(f'This analysis does not permit any conclusion about false-positives.')
