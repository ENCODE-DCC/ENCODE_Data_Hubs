__author__ = 'kelley'

import json
from general import convert_to_IHEC_format

VERSION='1.6'

def chip_seq_wrapper(assembly, taxon_id, target):
    url = 'https://www.encodeproject.org/search/?type=experiment&assay_term_name=ChIP-seq&target.name=%s-human' % target

    # Used to set is_main
    track_hierarchy = {'peak_calls': ['optimal idr thresholded peaks', 'conservative idr thresholded peaks',
                                'replicated peaks', 'peaks', 'hotspots'],
                           'signal': ['signal p-value', 'fold change over control', 'signal', 'raw signal']}

    def dataset_additions_f(experiment, json_object):

        #Set experiment_type
        json_object['experiment_attributes']['experiment_type'] = experiment['target']['label']

        return json_object

    return convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, dataset_additions_f)





if __name__ == "__main__":
    targets = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    for t in targets:
        data = chip_seq_wrapper(assembly='hg19', taxon_id=9606, target=t)
        with open('../output/%s_v%s.json' % (t, VERSION), 'w+') as outfile:
            json.dump(data, outfile, indent=4)
