__author__ = 'kelley'

import json
from general import convert_to_IHEC_format

def bisulfite_wrapper(version):
    url = 'https://www.encodeproject.org/search/?type=experiment&assay_term_name=whole-genome%20shotgun%20bisulfite%20sequencing'
    assembly = 'hg19'
    taxon_id = 9606

    # Used to set is_main
    track_hierarchy = {'methylation_profile': ['methylation state at CpG', 'methylation state at CHH']}

    def dataset_additions_f(experiment, json_object):

        #Set experiment_type
        json_object['experiment_attributes']['experiment_type'] = 'DNA Methylation'
        json_object['experiment_attributes']['assay_type'] = 'WGB-Seq'

        return json_object

    with open('../output/bisulfite_v%s.json' % version, 'w+') as outfile:
        json.dump(convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, dataset_additions_f), outfile, indent=4)

if __name__ == "__main__":

    ############################# Load bisulfite experiments #############################

    bisulfite_wrapper('1.6')