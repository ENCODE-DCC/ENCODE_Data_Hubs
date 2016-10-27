import requests
import json
from general import *
from rnaseq import *
from bisulfite import *
from chipseq import *



API_ENDPOINT = "https://www.encodeproject.org"
REFERENCE_EPIGENOME_COLLECTION = API_ENDPOINT + "/search/?type=ReferenceEpigenome&format=json&limit=all"

def run():
    hg19_dataset = dict()
    response = requests.get(API_ENDPOINT)
    
    # Iterate over Reference Epigenomes
    for reference_epigenome in response.json()['@graph']:
        reference_epigenome_object = requests.get(reference_epigenome['@id']+'?format=json').json()
        reference_registry_id = get_epirr_id(reference_epigenome_object)
        
        # Iterate over experiments in Reference Epigenomes
        for experiment_obj in reference_epigenome_object:
            addition = determine_addition(experiment_obj)
            if experiment_obj['assembly'][0] == 'hg19':
                hg19_dataset.update(create_datasets(experiment_obj, reference_registry_id, addition))
    
    with open('../output/hg19.json', 'w+') as outfile:
        json.dumps(hg19_dataset, outfile, indent=4)





if __name__ == "__main__":
    run()