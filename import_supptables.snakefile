# in order not to be dependent on running NGS pipelines, you can import supp tables.

rule all:
    input:  # Data from Supp. Tables 6 and 7 are not required for running any scripts
        'output/mrna_data_all.csv',
        'output/TableS2_RIP-seq.xlsx',
        'output/TableS5_Ribo-seq.xlsx',
        'output/TableS4_FullProteome.xlsx',
        'output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
        'output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
        # read count tables from GEO (TODO for this I need rules that download the files directly from the GEO website)
        'output/cluster_kos/read_counts.tsv',
        'output/sipool/read_counts.tsv',
        'output/rnai_read_counts.tsv',

rule mrna_data_all:
    input:
        'supp/TableS1_RNA-seq.xlsx'
    output:
        'output/mrna_data_all.csv',
    run:
        import pandas as pd
        df = pd.read_excel(input[0], skiprows=2, index_col=[0, 1], header=[0, 1])
        # source code compatibility
        df = df.rename(columns={'Ago2&1':'Ago12'})
        df.to_csv(output[0])

rule supp_table_s2:
    input:
        'supp/TableS2_RIP-seq.xlsx'
    output:
        'output/TableS2_RIP-seq.xlsx'
    shell: '''
        cp {input} {output}
    '''

rule ribo_seq:
    input:
        'supp/TableS5_Ribo-seq.xlsx'
    output:
        'output/TableS5_Ribo-seq.xlsx'
    run:
        import pandas as pd
        df = pd.read_excel(input[0], skiprows=2, index_col=[0, 1])
        df = df.rename(columns=lambda c: c.replace('AGO2&1', 'AGO12'))
        df.to_excel(output[0])

rule proteome:
    input:
        'supp/TableS4_FullProteome.xlsx'
    output:
        'output/TableS4_FullProteome.xlsx'
    run:
        import pandas as pd
        df = pd.read_excel(input[0], skiprows=2, index_col=0)
        df = df.rename(columns=lambda c: c.replace('AGO2&1', 'AGO12'))
        df.to_excel(output[0])

rule mirna_ko_quantseq:
    input:
        'supp/TableS6_miR-290-295KO_Quant-seq.xlsx'
    output:
        'output/cluster_kos/TableS6_miR-290-295KO_Quant-seq.xlsx',
    shell: '''
        cp {input} {output}
    '''
rule sipool_quantseq:
    input:
        'supp/TableS7_sipool-Tfap4-Quant-seq.xlsx'
    output:
        'output/sipool/TableS7_sipool-Tfap4-Quant-seq.xlsx',
    shell: '''
        cp {input} {output}
    '''

# TODO the following file needs to be downloaded from GEO (GSE181393)
rule quant_seq_read_counts:
    input:
        'supp/quant_seq_counts.tsv',
    output:
        cluster_kos='output/cluster_kos/read_counts.tsv',
        sipool='output/sipool/read_counts.tsv',
    run:
        import pandas as pd
        df = pd.read_csv(input[0], sep='\t', index_col=0)
        df.iloc[:, 0:4].to_csv(output['cluster_kos'], sep='\t')
        df.iloc[:, 4:8].to_csv(output['sipool'], sep='\t')


rule rnai_read_counts:
    input:
        'supp/TableS1_RNA-seq.xlsx'
    output:
        'output/rnai_read_counts.tsv',
    run:
        import pandas as pd
        df = pd.read_excel(input[0], skiprows=2, sheet_name='Raw read counts', index_col=[0,1])
        df = df.rename(columns=lambda c: c.replace('Ago2&1', 'Ago12'))
        df = df.droplevel(1)
        df.to_csv(output[0], sep='\t')
