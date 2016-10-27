import requests
import json
from general import *
from rnaseq import *
from bisulfite import *
from chipseq import *



API_ENDPOINT = "https://www.encodeproject.org"
REFERENCE_EPIGENOME_COLLECTION = API_ENDPOINT + "/search/?type=ReferenceEpigenome&format=json"

def collect_experiments(assembly, taxon_id):
    dataset = dict()
    count = 0
    response = requests.get(REFERENCE_EPIGENOME_COLLECTION)
    
    # Iterate over Reference Epigenomes
    for reference_epigenome in response.json()['@graph']:
        reference_epigenome_object = requests.get(API_ENDPOINT + reference_epigenome['@id']+'?format=json').json()
        reference_registry_id = get_epirr_id(reference_epigenome_object)

        # Iterate over experiments in Reference Epigenomes
        for experiment_obj in reference_epigenome_object['related_datasets']:
            if 'assembly' in experiment_obj and experiment_obj['assembly'] and experiment_obj['assembly'][0] == assembly:
                addition = determine_addition(experiment_obj)
                if addition:
                    dataset.update(create_datasets(experiment_obj, reference_registry_id, addition))
                    count = count + 1
                    print(count)

    data = convert_to_IHEC_format(dataset, assembly, taxon_id)
    return data