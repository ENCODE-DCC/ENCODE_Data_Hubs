__author__ = 'kelley'

import json
from general import convert_to_IHEC_format

VERSION='1.6'

def bisulfite_wrapper(assembly, taxon_id):
    url = 'https://www.encodeproject.org/search/?type=experiment&assay_term_name=whole-genome%20shotgun%20bisulfite%20sequencing'

    # Used to set is_main
    track_hierarchy = {'methylation_profile': ['methylation state at CpG', 'methylation state at CHH']}

    def dataset_additions_f(experiment, json_object):

        #Set experiment_type
        json_object['experiment_attributes']['experiment_type'] = 'DNA Methylation'
        json_object['experiment_attributes']['assay_type'] = 'WGB-Seq'

        return json_object

    return convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, dataset_additions_f)



if __name__ == "__main__":
    data = bisulfite_wrapper(assembly='hg19', taxon_id=9606)
    with open('../output/bisulfite_v%s.json' % VERSION, 'w+') as outfile:
        json.dump(data, outfile, indent=4)