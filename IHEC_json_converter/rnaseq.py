__author__ = 'kelley'

import json

VERSION='1.6'


RNASEQ_TRACK_HIEARCHY = {'signal_forward': ['plus strand signal of unique reads', 'plus strand signal of all reads',
                                              'plus strand signal', 'raw plus strand signal'],
                           'signal_reverse': ['minus strand signal of unique reads', 'minus strand signal of all reads',
                                              'minus strand signal', 'raw minus strand signal'],
                           'signal': ['signal of unique reads', 'signal of all reads', 'signal', 'raw signal', 'splice junctions'],
                           'contigs': ['contigs']}


def rna_seq_addition(experiment, json_object):

    size_range_to_experiment_type = {
        '>200': 'mRNA-seq',
        '<200': 'smRNA-seq'
    }

    #Set experiment_type
    size_range = None
    if 'size_range' in experiment['replicates'][0]['library']:
        size_range = experiment['replicates'][0]['library']['size_range']

    if size_range is None:
        print('Could not find size_range ' + experiment['accession'])
        json_object['experiment_attributes']['experiment_type'] = 'RNA-seq'
        json_object['experiment_attributes']['assay_type'] = 'RNA-seq'
    elif size_range not in size_range_to_experiment_type:
        print('Size range not found: ' + experiment['replicates'][0]['library']['size_range'])
        json_object['experiment_attributes']['experiment_type'] = 'RNA-seq'
        json_object['experiment_attributes']['assay_type'] = 'RNA-seq'
    else:
        json_object['experiment_attributes']['experiment_type'] = size_range_to_experiment_type[size_range]
        json_object['experiment_attributes']['assay_type'] = size_range_to_experiment_type[size_range]

    return json_object


# if __name__ == "__main__":
#     data = rna_seq_wrapper(assembly='hg19', taxon_id=9606)
#     with open('../output/RNAseq_v%s.json' % VERSION, 'w+') as outfile:
#         json.dump(data, outfile, indent=4)