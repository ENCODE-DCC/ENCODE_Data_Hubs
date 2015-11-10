__author__ = 'kelley'

import json
from general import convert_to_IHEC_format, set_main_track, signal_mapping

def rna_seq_wrapper(version):
    url = 'https://www.encodeproject.org/search/?type=experiment&assay_term_name=RNA-seq'
    assembly = 'hg19'
    taxon_id = 9606

    # Used to set is_main
    track_hierarchy = {'signal_forward': ['plus strand signal of unique reads', 'plus strand signal of all reads',
                                              'plus strand signal', 'raw plus strand signal'],
                           'signal_reverse': ['minus strand signal of unique reads', 'minus strand signal of all reads',
                                              'minus strand signal', 'raw minus strand signal'],
                           'signal': ['signal of unique reads', 'signal of all reads', 'signal', 'raw signal', 'splice junctions'],
                           'contigs': ['contigs']}

    size_range_to_experiment_type = {
        '>200': 'mRNA-seq',
        '<200': 'smRNA-seq'
    }

    def dataset_additions_f(experiment, json_object):

        #Set experiment_type
        size_range = None
        if 'size_range' in experiment['replicates'][0]['library']:
            size_range = experiment['replicates'][0]['library']['size_range']

        if size_range is None:
            print 'Could not find size_range ' + experiment['accession']
            json_object['experiment_attributes']['experiment_type'] = 'RNA-seq'
            json_object['experiment_attributes']['assay_type'] = 'RNA-seq'
        elif size_range not in size_range_to_experiment_type:
            print 'Size range not found: ' + experiment['replicates'][0]['library']['size_range']
            json_object['experiment_attributes']['experiment_type'] = 'RNA-seq'
            json_object['experiment_attributes']['assay_type'] = 'RNA-seq'
        else:
            json_object['experiment_attributes']['experiment_type'] = size_range_to_experiment_type[size_range]
            json_object['experiment_attributes']['assay_type'] = size_range_to_experiment_type[size_range]

        return json_object

    with open('../output/RNAseq_v%s.json' % version, 'w+') as outfile:
        json.dump(convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, dataset_additions_f), outfile, indent=4)

if __name__ == "__main__":

    ############################# Load RNA-seq experiments #############################

    rna_seq_wrapper('1.6')