__author__ = 'kelley'


VERSION='1.6'

# Used to set is_main
BISULFATE_TRACK_HIEARCHY = {'methylation_profile': ['methylation state at CpG', 'signal']}


def bisulfate_addition(experiment, json_object):
    #Set experiment_type
    json_object['experiment_attributes']['experiment_type'] = 'DNA Methylation'
    json_object['experiment_attributes']['assay_type'] = 'WGB-Seq'

    return json_object


# if __name__ == "__main__":
#     data = bisulfite_wrapper(assembly='hg19', taxon_id=9606)
#     with open('../output/bisulfite_v%s.json' % VERSION, 'w+') as outfile:
#         json.dump(data, outfile, indent=4)