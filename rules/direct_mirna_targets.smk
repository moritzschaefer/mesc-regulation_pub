from moritzsphd.util import remote_file
from snakemake.remote import HTTP
HTTP = HTTP.RemoteProvider()

rule gencode_gff_vM3:
    input:
        remote_file(config['gencode_gtf_vM3'], gunzip=True)
    output:
        'ref/gencode_vM3.db',
    conda:
        "../env/python.yaml"
    script:
        "../script/create_gffutils_db.py"

rule prepare_targetscan:
    input:
        conserved=HTTP.remote("http://www.targetscan.org/mmu_72/mmu_72_data_download/Conserved_Site_Context_Scores.txt.zip", static=True, keep_local=True),
        nonconserved=HTTP.remote("http://www.targetscan.org/mmu_72/mmu_72_data_download/Nonconserved_Site_Context_Scores.txt.zip", static=True, keep_local=True),
        annotation='ref/gencode_vM3.db',
    output: "output/targetscan_pairs.csv"
    threads: 4
    conda: "../env/python.yaml"
    script: "../script/combine_targetscan.py"


rule mesc_mirnas:
    input:
        mirna_cpms=HTTP.remote('https://fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/vis_viewer/data/mouse/mouse.srna.cpm.txt', static=True, keep_local=True),
        sample_info=HTTP.remote('https://fantom.gsc.riken.jp/5/suppl/De_Rie_et_al_2017/vis_viewer/data/mouse/mouse.srna.samples.tsv', static=True, keep_local=True)
    output:
        'output/mesc_mirnas.csv'
    log: 'log/mesc_mirnas.log'
    conda: '../env/python.yaml'
    script: '../script/mesc_mirnas.py'

rule ago2heap:
    'Use moritzsphd to get annotated list of heap-ago2 binding sites'
    output:
        ago2_heap_data='output/ago2_heap_mesc_mres.csv',
    threads: 1
    conda:
        '../env/python.yaml'
    script:
        '../script/prepare_ago2_heap_data.py'

rule combine_direct_mirna_interactions:
    '''
    returns a list of all potential pairs (targetscan) with additional fields:
    - heap peaks for this mirna
    - up in how many mutants
    - TS score
    - mirna expression
    '''
    input:
        mrna_data='output/mrna_data_{subset}.csv',
        ago2_heap_data='output/ago2_heap_mesc_mres.csv',
        targetscan_prediction="output/targetscan_pairs.csv",  # downloads.smk
        #  mirna_expression='output/TableS2_RIP-seq.xlsx',
        mirna_data='output/TableS2_RIP-seq.xlsx',
        annotation='ref/gencode.db', # TODO 'chr' should be added in prepare_targetscan
    output:
        main='output/mirnas/raw_mirna_interactions_{subset,(all|protein_coding)}.csv'
    params:
        mutants=config['full_effect_mutants'],
        combined_padj_threshold=config['combined_padj_threshold'],
        log2fc_threshold=config['log2fc_threshold'],
        min_mrna_expression=config['mrna_threshold']
    conda: '../env/python.yaml'
    script: '../script/direct_mirna_targets/combine_direct_mirna_interactions.py'

rule rank_direct_mirna_interactions:
    '''
    Rank, filter and group miRNA interactions
    '''
    input:
        raw_interactions='output/mirnas/raw_mirna_interactions_{subset}.csv',
        mrna_data='output/mrna_data_{subset}.csv',
    output:
        unfiltered='output/mirnas/unfiltered_interaction_ranking_{subset,(all|protein_coding)}.csv',
        filtered='output/mirnas/interaction_ranking_{subset}.csv',
        low_filtered='output/mirnas/low_interaction_ranking_{subset}.csv',
        grouped='output/mirnas/interaction_ranking_grouped_{subset}.csv', # includes filtering
        up_genes='output/up_genes_{subset}.txt',
        low_up_genes='output/low_up_genes_{subset}.txt',
        down_genes='output/down_genes_{subset}.txt'
    params:
        mutants=config['full_effect_mutants'],
        max_ts_score=config['ts_threshold'],
        min_mirna_expression=config['mirna_threshold'],
        min_num_up_genes=config['min_num_up_genes'],
        combined_padj_threshold=config['combined_padj_threshold'],
        log2fc_threshold=config['log2fc_threshold']
    conda: '../env/python.yaml'
    script: '../script/direct_mirna_targets/rank_direct_mirna_interactions.py'

rule plot_mirna_interaction_overlap:
    input:
        raw_interactions='output/mirnas/raw_mirna_interactions_{subset}.csv',
        up_genes="output/up_genes_{subset}.txt"
    output:
        overlap_plot='plot/mirna_interaction_overlap_{subset}.svg'
    conda: '../env/python.yaml'
    script: '../script/direct_mirna_targets/plot_mirna_interaction_overlap.py'

rule filter_esc_specific:
    'Split miRNA interactions based on miRNA-mES-specificity'
    input:
        mesc_mirnas='output/mesc_mirnas.csv',
        mirna_interaction_ranking='output/mirnas/interaction_ranking_{subset}.csv'
    output:
        mescspecific='output/mirnas/interaction_ranking_{subset}_mescspecific.csv',
        nonspecific='output/mirnas/interaction_ranking_{subset}_nonspecific.csv',
    conda: '../env/python.yaml'
    script: '../script/filter_esc_specific.py'

rule mirna_target_genes:
    input:
        'output/mirnas/interaction_ranking_grouped_all.csv'
    output:
        'output/mirnas/mirna_target_genes.txt'
    run:
        import pandas as pd
        with open(output[0], 'w') as f:
            f.write(','.join(pd.read_csv(input[0])['Geneid'].values))


GSEA_LIBRARIES = ['WikiPathways_2019_Mouse',
                  'GO_Biological_Process_2021',
                  'GO_Molecular_Function_2021',
                  'GO_Cellular_Component_2021',
                  'BioPlanet_2019',
                  'KEGG_2019_Mouse']
rule mirna_target_enrichment_analysis:
    input:
        libraries=expand('misc/{library}.gmt', library=GSEA_LIBRARIES),
        gene_sets=['output/mirnas/mirna_target_genes.txt',
                   'output/up_genes_all.txt',
                   'output/down_genes_all.txt'],
        enrichments=['output/direct_mirna_target_enrichment_all.csv',
                     'output/up_genes_enrichment_all.csv',
                     'output/down_genes_enrichment_all.csv'],
        human_mouse_trans='misc/human_mouse_trans.tsv'
    output:
        data='output/mirna_target_enrichment_analysis.xlsx',
        plot='plot/mirna_target_enrichment_analysis.svg'
    params:
        gene_set_labels=['miRNA target genes', 'Upregulated genes', 'Downregulated genes'],
        gsea_libraries=GSEA_LIBRARIES
    conda: '../env/python.yaml'
    script: '../script/mirna_target_enrichment_analysis.py'

rule analyze_direct_mirna_targets:
    '''
    Here we analyze direct miRNA targets.
    This mostly involves a GO analysis
    '''
    input:
        'output/mirnas/interaction_ranking_{subset}.csv',
    output:
        enrichment_table="output/direct_mirna_target_enrichment_{subset}.csv",
        enrichment_plot="plot/direct_mirna_target_enrichment_{subset}.svg",
    params:
        filter_padj=0.001,
        column='Gene name'
    conda: '../env/python.yaml'
    script: '../script/enrichment_analysis.py'

rule analyze_updown_genes:
    '''
    Here we analyze direct miRNA targets.
    This mostly involves a GO analysis
    '''
    input:
        'output/{direction}_genes_{subset}.txt',
    output:
        enrichment_table="output/{direction}_genes_enrichment_{subset}.csv",
        enrichment_plot="plot/{direction}_genes_enrichment_{subset}.svg",
    params:
        filter_padj=0.001,
        column=None
    conda: '../env/python.yaml'
    script: '../script/enrichment_analysis.py'

rule ribo_seq_validation_supptable:
    input:
        ribo_seq='/mnt/cclab_nas/groupdata/Articles Raw data/2021/Sch??fer et al. miRNA target prediction/copy_of_daniels/Tables/data_S4.xlsx',
    output:
        supp_table='output/TableS5_Ribo-seq.xlsx'
    run:
        import pandas as pd
        df_deg = pd.read_excel(input['ribo_seq'], sheet_name='RIBO_DESeq2', index_col=[0, 1], skiprows=2)
        df_deg = df_deg[[c for c in df_deg.columns if not c.startswith('cat_')]]
        df_counts = pd.read_excel(input['ribo_seq'], sheet_name='RIBO_spikeInNorm_counts', index_col=[0, 1], skiprows=2)
        writer = pd.ExcelWriter(output['supp_table'], engine='xlsxwriter')
        df_deg.to_excel(writer, sheet_name='Diff. expression analysis', startrow=1, index=True)
        writer.sheets['Diff. expression analysis'].write(0, 0, 'Supp. Table 4: Ribo-seq differential expression analysis of RNAi_KO mutants versus WT. Referenced by (Supp.) Figure 3')
        df_counts.to_excel(writer, sheet_name='spikeInNorm_counts', startrow=1, index=True)
        writer.sheets['spikeInNorm_counts'].write(0, 0, 'Drosophila-normalized Ribo-seq counts of RNAi_KO mutants and WT')
        writer.save()

rule ribo_seq_validation:
    input:
        ribo_seq='output/TableS5_Ribo-seq.xlsx',
        interaction_ranking='output/mirnas/interaction_ranking_grouped_{group}.csv'
    output:
        plot='plot/ribo_seq_validation_{group}.svg',
    params:
        mutants=config['full_effect_mutants'],
        sample_colors=config['sample_colors']
    log: 'log/ribo_seq_validation_{group}.txt'
    conda: '../env/python.yaml'
    script: '../script/ribo_seq_validation.py'

rule ms_validation_supptable:
    input:
        ms='/mnt/cclab_nas/groupdata/Articles Raw data/2021/Sch??fer et al. miRNA target prediction/copy_of_daniels/Tables/data_S3.xlsx',
    output:
        supp_table='output/TableS4_FullProteome.xlsx'
    run:
        import pandas as pd
        df_deg = pd.read_excel(input['ms'], sheet_name='MS_DE', skiprows=2)
        df_counts = pd.read_excel(input['ms'], sheet_name='MS_log2Quant', skiprows=2)
        writer = pd.ExcelWriter(output['supp_table'], engine='xlsxwriter')
        df_deg.to_excel(writer, sheet_name='Diff. expression analysis', startrow=1, index=False)
        writer.sheets['Diff. expression analysis'].write(0, 0, 'Supp. Table 5: differential expression analysis for protein abundances of RNAi_KO mutants versus WT. Referenced by (Supp.) Figure 3')
        df_counts.to_excel(writer, sheet_name='log2_quantification', startrow=1, index=False)
        writer.sheets['log2_quantification'].write(0, 0, 'log2-quantification of protein abundance for RNAi_KO mutants and WT. samples with b extension are technical replicates, proteins with a prot.FDR > 0.02 were filtered.')
        writer.save()

rule ms_validation:
    input:
        ms='output/TableS4_FullProteome.xlsx',
        interaction_ranking='output/mirnas/interaction_ranking_grouped_all.csv'
    output:
        plot='plot/ms_validation_all.svg',
    params:
        mutants=config['full_effect_mutants'],
        sample_colors=config['sample_colors']
    log: 'log/ms_validation_all.txt'
    conda: '../env/python.yaml'
    script: '../script/ms_validation.py'

rule target_cluster_distribution:
    input:
        mirna_data='output/TableS2_RIP-seq.xlsx',
        interaction_ranking='output/mirnas/interaction_ranking_{group}.csv'
    output:
        barplot='plot/target_cluster_distribution_barplot_{group}.svg',
        pieplot='plot/target_cluster_distribution_pieplot_{group}.svg'
    conda: '../env/python.yaml'
    script: '../script/target_cluster_distribution.py'

rule cooperative_mirna_binding:
    input:
        interaction_ranking='output/mirnas/interaction_ranking_all.csv',
        unfiltered_interactions='output/mirnas/unfiltered_interaction_ranking_all.csv',
        ribo_data='output/TableS5_Ribo-seq.xlsx',
        quantseq_data='output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        ago2_heap_peaks=remote_file('https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE139345&format=file&file=GSE139345_mESC_peaks.csv.gz', gunzip=True)
    output:
        heap_for_num_mres='plot/heap_for_num_mres.svg',
        interaction_count_ribo='plot/interaction_count_ribo.svg',
        interaction_count_mir290='plot/interaction_count_mir290.svg'
    params:
        mre_count_bins=config['mre_count_bins'],
        mutants=config['full_effect_mutants'],
        ts_mre_count_bins=config['ts_mre_count_bins'],
        mirna_threshold=config['mirna_threshold'],
        sample_colors=config['sample_colors']
    conda: '../env/python.yaml'
    script:
        '../script/cooperative_mirna_binding.py'

rule database_overlap:
    input:
        interaction_ranking='output/mirnas/interaction_ranking_all.csv',
        mtb_mmu=HTTP.remote('https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/cache/download/8.0/mmu_MTI.xls', keep_local=True, static=True),
        mirna_data='output/TableS2_RIP-seq.xlsx',
    output:
        plot='plot/mirtarbase_overlap.svg',
        data='output/mirtarbase_overlap.xlsx'
    params:
        mirna_threshold=config['mirna_threshold']
    conda: '../env/python.yaml'
    script:
        '../script/database_overlap.py'

rule interaction_stats:
    input:
        interaction_ranking='output/mirnas/interaction_ranking_all.csv',
        mrna_data='output/mrna_data_all.csv',
    output:
        interaction_histo='plot/interaction_count_histo.svg',
        gene_expression='plot/interaction_count_gene_expr.svg',
        mrna_sections='plot/interaction_count_mrna_sections.svg',
        mre_type='plot/interaction_count_mre_type.svg'
    params:
        mre_count_bins=config['mre_count_bins'],
    conda: '../env/python.yaml'
    script:
        '../script/interaction_stats.py'

