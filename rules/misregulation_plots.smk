'''
Mostly figure 1 of my paper
'''
rule ma_plots:
    input: 'output/mrna_data_{subset}.csv'
    output: 'plot/ma_{subset}.{ext}'
    params:
        mutants=config['full_effect_mutants']
    conda: '../env/python.yaml'
    script: '../script/ma_plots.py'

rule misregulation_overlap_plots:
    input: 'output/mrna_data_{subset}.csv'
    output: 'plot/overlaps_{subset}_{direction}.svg',
    params:
        mutants=config['full_effect_mutants']
    conda: '../env/eulerr.yaml'  # doesn't work in nixos+conda
    script: '../script/misregulation_overlap_plots.R'

rule misregulation_pairplots:
    input: 'output/mrna_data_{subset}.csv'
    output: 'plot/misregulation_pairplot{filtered,(_filtered)?}_{subset}.svg'
    params:
        mutants=config['full_effect_mutants']
    conda: '../env/python.yaml'
    script: '../script/misregulation_pairplots.py'

rule misreg_parameter_selection:
    'Show which parameters lead to what numbers of misregulation'
    input: 'output/mrna_data_{subset}.csv'
    output:

# PCA plots
rule mrna_pca_plot:
    input:
        expression='output/rnai_read_counts.tsv',
    params:
        samples=FULL_EFFECT_WT_SAMPLES,
        sample_colors=config['sample_colors']
    conda: '../env/python.yaml'
    output: 'plot/mrna_pca.{ext}'
    script: '../script/mrna_pca_plot.py'

rule mirna_pca_plot:
    input:
        mirna_data='output/TableS2_sRNA-seq.xlsx'
    params:
        samples=FULL_EFFECT_WT_SAMPLES,
        sample_colors=config['sample_colors']
    conda: '../env/python.yaml'
    output: 'plot/mirna_pca.{ext}'
    script: '../script/mirna_pca_plot.py'

rule mirna_cdf_plot:
    input:
        'output/mirna_data.csv'
    params:
        mutants=config['full_effect_mutants'],
        sample_colors=config['sample_colors']
    conda: '../env/python.yaml'
    output: 'plot/mirna_cdf.{ext}'
    script: '../script/mirna_cdf_plot.py'


# Correlation plot
rule transcript_correlation_plot:
    input:
        mrna_data='output/mrna_data_all.csv',
        expression='output/rnai_read_counts.tsv',
    params:
        samples=FULL_EFFECT_WT_SAMPLES,
        mutants=config['full_effect_mutants']
    conda: '../env/python.yaml'
    output: 'plot/transcript_correlation.{ext}'
    script: '../script/transcript_correlation_plot.py'

rule common_genes_matrix:
    input:
        mrna_data='output/mrna_data_all.csv',
    params:
        mutants=config['full_effect_mutants']
    conda: '../env/python.yaml'
    output: 'plot/gene_{direction}regulation.{ext}'
    script: '../script/gene_upregulation_plot.py'
