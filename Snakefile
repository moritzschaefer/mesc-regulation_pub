import itertools
from moritzsphd.util import remote_file
##### load config and sample sheets #####
configfile: "config.yaml"

ALL_MUTANTS = config['full_effect_mutants'] + config['other_mutants']
FULL_EFFECT_WT_SAMPLES = [v for sub in [(f'{s}_1', f'{s}_2') for s in config['full_effect_mutants'] + ['WT']] for v in sub]
SAMPLES = [v for sub in [(f'{s}_1', f'{s}_2') for s in ALL_MUTANTS + ['WT']] for v in sub]

include: "rules/tpm.smk"
include: 'rules/direct_mirna_targets.smk'
include: 'rules/misregulation_plots.smk'
include: 'rules/tf_annotation.smk'
include: 'rules/cluster_kos.smk'
include: 'rules/sipool.smk'
include: 'rules/figure_plots.smk'
include: 'rules/clip_plots.smk'


include: 'rules/figures.smk'  # in this file, figure-files are defined on a per-figure basis
rule all:
    input:
        rules.figure1.input,
        rules.figure2.input,
        rules.figure3.input,
        rules.figure4.input


# TODO maybe delete this rule? the supp. tables were generated manually in any case.
# The order of the supp tables might have changed in the final publication. Also, some sheets might have been added manually
rule supp_tables:
    input:
        # 'output/TableS1_RNA-seq.xlsx',  # generated manually
        # 'output/TableS2_RIP-seq.xlsx',  # generated manually
        'output/TableS3_Integrative_analysis.xlsx',
        'output/TableS5_Ribo-seq.xlsx',
        'output/TableS4_FullProteome.xlsx',
        'output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        'output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
        # read count tables from GEO
        'output/cluster_kos/read_counts.tsv',
        'output/sipool/read_counts.tsv',
        'output/rnai_read_counts.tsv',


## RNA-seq
rule rnai_read_count_matrix:
    input:
        remote_file(config['expression_path'])
    output:
        read_count_matrix='output/rnai_read_counts.tsv'
    params:
        samples=FULL_EFFECT_WT_SAMPLES
    run:
        import pandas as pd

        df = pd.read_csv(input[0], sep='\t', index_col=0)
        df[params['samples']].to_csv(output[0], sep='\t')

# Combine the RNA-seq data. commented, because this data is imported from supp. table
# rule combine_mrna_data:
#     input:
#         tpm_expression='output/all.tpm.tsv',
#         diffexp=[remote_file(config['mrna_diffexp_path'].format(sample=mutant)) for mutant in config['full_effect_mutants']]
#     params:
#         mutants=config['full_effect_mutants']
#     conda: 'env/python.yaml'
#     output:
#         csv='output/mrna_data_all.csv',
#         supp_table='output/TableS1_RNA-seq_rna_seq.xlsx',
#     script: 'script/combine_mrna_data.py'

rule filter_protein_coding:
    '''
    Extract protein coding genes
    also log a comparison of
    how many NC-genes were affected as compared to coding genes
    '''
    input:
        mrna_data='output/mrna_data_all.csv',
        annotation='ref/gencode.db'
    log: 'log/filter_protein_coding.log'
    conda: 'env/python.yaml'
    output: 'output/mrna_data_protein_coding.csv'
    script: 'script/filter_protein_coding.py'

rule gencode_gff:
    input:
        remote_file(config['gencode_gtf'], gunzip=True)
    output:
        "ref/gencode.db"
    conda:
        "env/python.yaml"
    script:
        "script/create_gffutils_db.py"
