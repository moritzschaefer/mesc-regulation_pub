'''Note: Some of the rules here might require additional processing steps which are not covered by this pipeline but described in the paper's Methods section'''
from snakemake.remote import HTTP
HTTP = HTTP.RemoteProvider()

rule mesc_tfap4_genome_track:
    input:
        tracks='misc/mesc_tfap4_tracks.ini',
        clip_track='misc/bosson_tagged_Ago2_iCLIP.bw.mm10',
        heap_track='output/track_plots/heap_peaks.bed',
        expressed_mres='output/track_plots/expressed_mres.bed'
    output:
        'plot/mesc_tfap4_track.pdf',
    shell: '''
        pyGenomeTracks --tracks {input.tracks} --region 16:4544545-4545861 --outFileName {output} --plotWidth 20 --trackLabelFraction '0.0'
    '''

rule download_annotation:
    input:
        HTTP.remote('https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz', static=True)
    output:
        'misc/human_anno.gtf'
    shell: '''
    gunzip < {input} > {output}'''

rule hesc_tfap4_genome_track:
    '''This is the hESC CLIP in which I look for Tfap4<->302c-3p binding
    all associated data is here: /home/moritz/wiki/roam/data/6f/be6c20-456c-4aff-b630-44b40dd3eba6

    302c-3p 8mer binding site in human Tfap4 (https://www.ensembl.org/Homo_sapiens/Transcript/Sequence_cDNA?db=core;g=ENSG00000090447;r=16:4257182-4258054;t=ENST00000204517) is at position ~442->450 (counted from the back...) (https://www.targetscan.org/cgi-bin/vert_80/view_gene.cgi?gs=TFAP4&taxid=9606&members=miR-302c-3p.2/520-3p&showcnc=1&shownc=1)

    mRNA-seq  --local  -i . -o out/ -c out/mRNA-seq.config.yaml GRCh38_gencode_release29
    '''

    input:
        tracks='misc/hesc_tfap4_tracks.ini',
        mre_track='misc/human_tfap4_mir302c_mre.bed',
        human_anno='misc/human_anno.gtf'
    output:
        'plot/hesc_clip_tfap4_track.pdf',
    shell: '''
    pyGenomeTracks --tracks {input.tracks} --region 16:4256965-4258280 --outFileName {output} --plotWidth 15 --trackLabelFraction '0.0'
    '''

rule hesc_interaction_conservation:
    '''
    The discrepancy between the control and the positive group clearly indicates a biological relevance of human-conserved interactions. Focusing on the set of conserved interactions with both hESC miRNA expression and associated hESC CLIP-seq reads (src_python[:session py]{len(df.loc[(df['human_mirna_family_cpm'] > 10) & (df['human_peak_size'] > 0)])} {{{results(=858=)}}} interactions), we also see a weak, but clear correlation between hESC miRNA expression and CLIP-seq read counts.
    '''
    input:
        hesc_clip_bw='misc/hesc_clip_SRR359787.coverage.bw',
        interaction_ranking='output/mirnas/interaction_ranking_all.csv',
        orthologs='misc/human_mouse_orthologs.tsv',
        hesc_mirna_expr='misc/hinton14_tables2.xlsx'
    output:
        mre_count_with_clip_reads='plot/clip_plots/mre_count_with_clip_reads.svg',
        mirna_expr_clip_conservation='plot/clip_plots/hesc_clip_conservation.svg',
        mirna_expr_clip_correlation='plot/clip_plots/hesc_interaction_conservation.svg',
        data='output/hesc_conservation_plot.xlsx'
    script:
        '../script/clip_plots/hesc_interaction_conservation.py'
