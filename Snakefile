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

# The order of the supp tables might have changed in the final publication. Also, some sheets might have been added manually
rule supp_tables:
    input:
        'output/TableS1_RNA-seq.xlsx',
        'output/TableS2_sRNA-seq.xlsx',
        'output/TableS3_Integrative_analysis.xlsx',
        'output/TableS5_Ribo-seq.xlsx',
        'output/TableS4_FullProteome.xlsx',
        'output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        'output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
        # read count tables from GEO
        'output/cluster_kos/read_counts.tsv',
        'output/sipool/read_counts.tsv',
        'output/rnai_read_counts.tsv',

# TODO fix all of them
rule figure1:
    input:
        'plot/misregulation_pairplot_all.svg',
        'plot/mirna_cdf.svg',
        'plot/ma_all.svg',
        'plot/mirna_pca.svg',
        'plot/mrna_pca.svg',
        'plot/transcript_correlation.svg',
        'plot/gene_upregulation.svg',
        'plot/gene_downregulation.svg'

rule figure2:
    input:
        expand('plot/target_data/{plot}_{gene}.pdf', plot=['track', 'target_log2fc'], gene=['Tfap4', 'Axin2', 'Rps26', 'Ctcf', 'Apoe',]), # TODO e.g.
        'plot/feature_combination_count_horizontal.svg'

rule figure3:
    input:
        'plot/ribo_seq_validation_all.svg',
        'plot/ms_validation_all.svg',
        'plot/target_cluster_distribution_barplot_all.svg',
        'plot/cluster_kos/ma_miR-290-295_all.svg',
        'plot/cluster_kos/target_cdf_miR-290-295.svg',

rule figure4:
    input:
        'plot/sipool/target_rescue.svg',
        'plot/mir_290_regulated_tfs.svg'

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

# Combine the RNA-seq
rule combine_mrna_data:
    input:
        tpm_expression='output/all.tpm.tsv',
        diffexp=[remote_file(config['mrna_diffexp_path'].format(sample=mutant)) for mutant in config['full_effect_mutants']]
    params:
        mutants=config['full_effect_mutants']
    conda: 'env/python.yaml'
    output:
        csv='output/mrna_data_all.csv',
        supp_table='output/TableS1_RNA-seq_rna_seq.xlsx',
    script: 'script/combine_mrna_data.py'

# relies on moritzsphd.data
rule combine_mirna_data:
    conda: 'env/python.yaml'
    output:
        csv='output/mirna_data.csv',
        supp_table='output/TableS2_sRNA-seq.xlsx'
    params:
        mutants=config['full_effect_mutants']
    script: 'script/combine_mirna_data.py'

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
