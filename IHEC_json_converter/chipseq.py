__author__ = 'kelley'

import json
from general import convert_to_IHEC_format

def chip_seq_wrapper(target, version):
    url = 'https://www.encodeproject.org/search/?type=experiment&assay_term_name=ChIP-seq&target.name=%s-human' % target
    assembly = 'hg19'
    taxon_id = 9606

    # Used to set is_main
    track_hierarchy = {'peak_calls': ['optimal idr thresholded peaks', 'conservative idr thresholded peaks',
                                'replicated peaks', 'peaks', 'hotspots'],
                           'signal': ['signal p-value', 'fold change over control', 'signal', 'raw signal']}

    def dataset_additions_f(experiment, json_object):

        #Set experiment_type
        json_object['experiment_attributes']['experiment_type'] = experiment['target']['label']

        return json_object

    with open('../output/%s_v%s.json' % (target, version), 'w+') as outfile:
        json.dump(convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, dataset_additions_f), outfile, indent=4)

if __name__ == "__main__":

    ############################# Load Chip-seq experiments #############################

    chip_seq_wrapper('H3K27ac', '1.6')
    chip_seq_wrapper('H3K27me3', '1.6')
    chip_seq_wrapper('H3K36me3', '1.6')
    chip_seq_wrapper('H3K4me1', '1.6')
    chip_seq_wrapper('H3K4me3', '1.6')
    chip_seq_wrapper('H3K9me3', '1.6')