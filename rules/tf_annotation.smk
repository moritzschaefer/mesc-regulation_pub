from snakemake.remote import HTTP
HTTP = HTTP.RemoteProvider()

rule define_tfs:
    '''
    - download TF data set from mTFkb
    - filter for expressed TFs (in WT or any of the mutants)
    - annotate TF expression (expressed only in ESC or also in others)
    '''
    input:
        tf_data=HTTP.remote('https://sunlab.cpy.cuhk.edu.hk/mTFkb/download/mTFkb.tar.gz', keep_local=True, static=True),
        mrna_data='output/mrna_data_protein_coding.csv'
    params:
        add_tfs=config['add_tfs'],
        remove_tfs=config['remove_tfs']
    output: 'output/tf_annotation.csv'
    conda: '../env/python.yaml'
    script: '../script/define_tfs.py'

rule annotate_mirna_ranking_tfs:
    '''
    Extract TFs from gene sets
    I have several lists, so in the output, I have TFs in indices and as columns I have information about different gene sets (e.g. identified in mirna_interaction_ranking_all.csv)

    '''
    input:
        tf_annotation='output/tf_annotation.csv',
        data='output/mirnas/interaction_ranking_{subset}.csv'
    output:
        data='output/mirnas/interaction_ranking_{subset}_istf.csv'
    run:
        import pandas as pd

        tfs = pd.read_csv(input['tf_annotation'], index_col=0)
        df = pd.read_csv(input['data'], index_col=0)
        df.insert(0, 'is_tf', df.index.map(lambda ensid: ensid in tfs.index))

        df.to_csv(output['data'])

rule annotate_ml_prediction_tfs:
    input:
        tf_annotation='output/tf_annotation.csv',
        data='output/mir290_ml/trained/linear_targets_outer_nonweighted.csv' # mirnas/interaction_ranking_{subset}
    output:
        data='output/mir290_ml/trained/linear_targets_outer_nonweighted_istf.csv'
    run:
        import pandas as pd

        tfs = pd.read_csv(input['tf_annotation'], index_col=1)
        df = pd.read_csv(input['data'], index_col=0, header=None)
        df.insert(0, 'is_tf', df.index.map(lambda ensid: ensid in tfs.index))
        print(df['is_tf'].sum())

        df.to_csv(output['data'])
