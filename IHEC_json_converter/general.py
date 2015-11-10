f__author__ = 'kelley'


import sys
import urlparse
import requests
import json
import os

'''
This file contains methods that convert ENCODE data (from its webservice urls) into JSON that IHEC then loads.
I tried to separate general code from assay-specific code in order to keep things fairly clean. This file contains
general-purpose methods. The assay-specific functions are in several other files in this directory. This code is meant
to be run from the main methods of the assay-specific files.
'''

signal_mapping = {
                  #Chip Seq
                  'raw signal': 'signal',
                  'signal': 'signal',
                  'signal p-value': 'signal',
                  'fold change over control': 'signal',
                  'optimal idr thresholded peaks': 'peak_calls',
                  'conservative idr thresholded peaks': 'peak_calls',
                  'replicated peaks': 'peak_calls',
                  'peaks': 'peak_calls',
                  'hotspots': 'peak_calls',

                  #RNA Seq
                  'raw plus strand signal': 'signal_forward',
                  'raw minus strand signal': 'signal_reverse',
                  'plus strand signal': 'signal_forward',
                  'minus strand signal': 'signal_reverse',
                  'contigs': 'contigs',
                  'plus strand signal of unique reads': 'signal_forward',
                  'minus strand signal of unique reads': 'signal_reverse',
                  'plus strand signal of all reads': 'signal_forward',
                  'minus strand signal of all reads': 'signal_reverse',
                  'signal of all reads': 'signal',
                  'signal of unique reads': 'signal',
                  'splice junctions': 'signal',

                  #Methylation
                  'methylation state at CpG': 'methylation_profile',
                  'methylation state at CHH': 'methylation_profile',
}

biomaterial_mapping = {
    'immortalized cell line': 'Cell Line',
    'tissue': 'Primary Tissue',
    'primary cell': 'Primary Cell',
    'stem cell': 'Cell Line',
    'in vitro differentiated cells': 'Primary Cell',
    'induced pluripotent stem cell line': 'Cell Line'
}

############################# General Code #############################
# This code is used to generate IHEC json objects from any type of ENCODE assay. Anything assay-specific belongs
# in other files.

def convert_to_IHEC_format(url, assembly, taxon_id, track_hierarchy, assay_specific_additions, limit='all'):
    #Prepare the url - make sure it has the right parameters.
    url = prep_query_url(url, assembly, taxon_id, limit=limit)

    #Get json response
    json_response = get_json_from_encode(url)

    if len(json_response) == 0:
        raise Exception('This url returns no data.')

    #Create hub description
    hub_description = create_hub_description(assembly, taxon_id)

    #Create all of the datasets (one for each experiment, sample_id pair)
    print '%d experiments returned from this query.' % len(json_response)

    datasets = dict()
    for entry in json_response:
        datasets.update(create_datasets(entry, assay_specific_additions))

    print '%d IHEC datasets created.' % len(datasets)

    # Merge datasets. We can merge two datasets if either
    # 1. they have the same experiment_id and their sample_attributes are identical
    # 2. they have the same experiment_type and the same biosample_id and they are from ENCODE3

    # Merge 1
    experiment_id_to_datasets = dict()
    for dataset in datasets.values():
        experiment_id = dataset['accession']
        if experiment_id not in experiment_id_to_datasets:
            experiment_id_to_datasets[experiment_id] = []
        experiment_id_to_datasets[experiment_id].append(dataset)

    datasets = dict()
    for experiment_id, these_datasets in experiment_id_to_datasets.iteritems():
        if is_match_sample_attributes([x['sample_attributes'] for x in these_datasets]):
            print 'Match found: ' + experiment_id

            dataset = these_datasets[0]
            dataset['browser'] = merge_tracks(these_datasets)
            dataset['sample_attributes']['sample_id'] = '-'.join(sorted([x['sample_attributes']['sample_id'] for x in these_datasets]))
            datasets[experiment_id] = dataset
        else:
            for dataset in these_datasets:
                datasets[experiment_id + '-' + dataset['sample_attributes']['sample_id']] = dataset

    # Merge 2
    group_key_to_datasets = dict()
    nogroup_datasets = dict()
    for accession, dataset in datasets.iteritems():
        if dataset['award'] == 'ENCODE3':
            key = (dataset['experiment_attributes']['experiment_type'], dataset['sample_attributes']['sample_id'])
            if key not in group_key_to_datasets:
                group_key_to_datasets[key] = []
            group_key_to_datasets[key].append((accession, dataset))
        else:
            nogroup_datasets[accession] = dataset

    datasets = dict(nogroup_datasets)
    for key, these_datasets in group_key_to_datasets.iteritems():
         if len(these_datasets) > 1:
            print 'Collapsed: ', key, len(these_datasets)

            dataset = these_datasets[0][1]
            dataset['browser'] = collapse_tracks([x[1] for x in these_datasets])
            datasets['-'.join([x[0] for x in these_datasets])] = dataset
         else:
             datasets[these_datasets[0][0]] = these_datasets[0][1]

    #Set is_main on tracks
    for accession, dataset in datasets.iteritems():
        for track_type in set(signal_mapping.values()):
            if track_type in dataset['browser']:
                set_main_track(dataset['browser'][track_type], track_hierarchy[track_type], accession, track_type)


    # Remove extra entries:
    for accession, dataset in datasets.iteritems():
        del dataset['award']
        del dataset['accession']
        del dataset['sample_attributes']['replicate']

    print '%d IHEC datasets created after merge.' % len(datasets)

    return {
        'hub_description': hub_description,
        'datasets': datasets
    }


def merge_tracks(datasets):
    # We want to merge these datasets together, so we need to merge tracks.
    new_tracks = {}
    for dataset in datasets:
        for track_type, tracks in dataset['browser'].iteritems():
            if track_type not in new_tracks:
                new_tracks[track_type] = []

            for track in tracks:
                track['subtype'] = track['subtype'] + ' (rep' + str(dataset['sample_attributes']['replicate']) + ')'
            new_tracks[track_type].extend(tracks)

    #If a track has been assigned to multiple replicates, we might end up seeing it twice - check for that
    cleaned_tracks = dict()
    for track_type, tracks in new_tracks.iteritems():
        url_to_track = dict()
        for track in tracks:
            if track['big_data_url'] not in url_to_track:
                url_to_track[track['big_data_url']] = []
            url_to_track[track['big_data_url']].append(track)

        for url, url_tracks in url_to_track.iteritems():
            if len(url_tracks) > 1:
                track = url_tracks[0]
                track['subtype'] = track['subtype'][0:track['subtype'].index('(rep')-1]
                url_to_track[url] = [track]
        cleaned_tracks[track_type] = [x[0] for x in url_to_track.values()]
    return cleaned_tracks

def collapse_tracks(datasets):
    # We want to collapse these datasets together, so we need to collapse tracks.
    new_tracks = {}
    for dataset in datasets:
        for track_type, tracks in dataset['browser'].iteritems():
            if track_type not in new_tracks:
                new_tracks[track_type] = []
            new_tracks[track_type].extend(tracks)
    return new_tracks

def create_hub_description(assembly, taxon_id):
    return {
        'taxon_id': taxon_id,
        'assembly': assembly,
        'publishing_group': 'ENCODE',
        'email': 'encode-help@lists.stanford.edu'
    }


def create_datasets(experiment, assay_specific_additions):
    print experiment['accession']
    experiment = get_json_from_encode('https://www.encodeproject.org/experiments/%s/?format=json' % experiment['accession'])

    #Create sample_attributes
    sample_attributes = [create_sample_attribute(replicate) for replicate in experiment['replicates']]

    #Create experiment_attributes
    experiment_attributes = create_experiment_attributes(experiment)

    #Create tracks
    datasets = dict()

    for sample_attribute in sample_attributes:
        replicate_num = sample_attribute['replicate']

        #Grab protocol document
        protocol_document = None
        try:
            protocol_document = 'documents/%s/%s' % (experiment['replicates'][0]['library']['documents'][0]['uuid'], experiment['replicates'][0]['library']['documents'][0]['attachment']['href'])
        except:
            'No protocol document for %s %d' % (experiment['accession'], replicate_num)

        tracks = filter(None, [create_track(file, protocol_document, replicate_num) for file in experiment['files']])

        all_tracks = {}
        for track_type in set(signal_mapping.values()):
            these_tracks = [track for key, track in tracks if key == track_type]
            if len(these_tracks) > 0:
                all_tracks[track_type] = these_tracks

        json_object = {
            'sample_attributes': sample_attribute,
            'experiment_attributes': experiment_attributes,
            'browser': all_tracks,

            # This is put in here temporarily in order to find datasets that can be merged into a single dataset
            # since they contain the same biosamples. We can only do this for ENCODE3 experiments, since the loading of
            # ENCODE3 chipseq experiments doesn't allow for this type of merging.
            'award': experiment['award']['rfa'],
            'accession': experiment['accession']
        }

        # Add any assay-specific stuff
        assay_specific_additions(experiment, json_object)
        datasets['%s-%s' % (experiment['accession'], json_object['sample_attributes']['sample_id'])] = json_object

    return datasets

def set_main_track(tracks, track_hierarchy, exp_accession, track_type):
    main_track = None
    main_track_hier = None
    for track in tracks:
        #Looking for main track
        subtype = track['subtype']
        replicate_offset = 0
        if '(rep' in subtype:
            replicate_offset = .1*int(subtype[subtype.index('(rep')+4:-1])
            subtype = subtype[0:subtype.index('(rep')-1]

        hierarchy = track_hierarchy.index(subtype)+replicate_offset
        if main_track_hier is None or hierarchy < main_track_hier:
            main_track = track
            main_track_hier = hierarchy

    #Set main track
    if main_track is None:
        print 'Dataset %s has no %s tracks.' % (exp_accession, track_type)
    else:
        main_track['is_main'] = True
    return tracks


def create_track(file, protocol_document_href, replicate_num):
    if file['file_format'] == 'bigWig' or file['file_format'] == 'bigBed':

        if file['output_type'] in signal_mapping:

            if replicate_num in file['biological_replicates'] or len(file['biological_replicates']) == 0:
                return signal_mapping[file['output_type']], {
                        'big_data_url': 'https://www.encodeproject.org%s@@download/%s.%s?proxy=true' % (file['@id'], file['accession'], file['file_format']),
                        'description_url': '' if protocol_document_href is None else 'https://www.encodeproject.org/%s?proxy=true' % protocol_document_href,
                        "md5sum": file['md5sum'],
                        "is_main": False,
                        "subtype": file['output_type']
                    }
        else:
            print 'Output type not found: ' + file['output_type']


def is_match_sample_attributes(sample_attributes):
    #Check that information matches among sample_attributes
    for i, x in enumerate(sample_attributes):
        for j, y in enumerate(sample_attributes[i+1:]):
            mismatching_entries = set(x.items()) - set(y.items())
            #It's ok if sample_ids and replicate numbers don't match.
            if len(mismatching_entries) > 2:
                return False
    return True


def create_sample_attribute(replicate):
    biosample = replicate['library']['biosample']

    sample_ontology_uri = None
    if biosample['biosample_term_id'].startswith('CL') or biosample['biosample_term_id'].startswith('UBERON'):
        sample_ontology_uri = 'http://purl.obolibrary.org/obo/%s' % biosample['biosample_term_id'].replace(':', '_')
    elif biosample['biosample_term_id'].startswith('EFO'):
        sample_ontology_uri = 'http://www.ebi.ac.uk/efo/%s' % biosample['biosample_term_id'].replace(':', '_')
    elif biosample['biosample_term_id'].startswith('NTR'):
        sample_ontology_uri = ''
    else:
        print 'Could not find url for sample ontology %s' % biosample['biosample_term_id']

    sample_attribute = {
        'sample_id' : biosample['accession'],
        'sample_ontology_uri' : sample_ontology_uri,
        'molecule' : replicate['library']['nucleic_acid_term_name'],
        'molecule_ontology_uri' : 'http://www.sequenceontology.org/miso/current_svn/term/%s' % replicate['library']['nucleic_acid_term_id'],
        'disease' : biosample.get('health_status', ''),
        'biomaterial_type' : biomaterial_mapping[biosample['biosample_type']],
        'disease_ontology_uri': '',
        'replicate': replicate['biological_replicate_number']
    }

    if 'donor' in biosample:
        add_SA_with_donor(sample_attribute, biosample['donor'])

    add_SA_cell_line(sample_attribute, biosample)
    add_SA_primary_cell(sample_attribute, biosample)
    #add_SA_primary_cell_culture(sample_attribute, biosample)
    add_SA_primary_tissue(sample_attribute, biosample)
    return sample_attribute


def add_SA_with_donor(sample_attribute, donor):
    sample_attribute['donor_id'] = donor['accession']
    sample_attribute['donor_age'] = donor.get('age', '')
    sample_attribute['donor_age_units'] = donor.get('age_units', '')
    sample_attribute['donor_life_stage'] = donor.get('life_stage', '')
    sample_attribute['donor_health_status'] = donor.get('health_status', '')
    sample_attribute['donor_sex'] = donor.get('sex', '').title()
    sample_attribute['donor_ethnicity'] = donor.get('ethnicity', '').title()

def add_SA_cell_line(sample_attribute, biosample):
    if sample_attribute['biomaterial_type'] == 'Cell Line':
        sample_attribute['line'] = biosample['biosample_term_name']
        sample_attribute['differenciate_stage'] = biosample['biosample_term_name']
        sample_attribute['sex'] = biosample['sex'].title()

        #Get lineage
        if biosample['biosample_type'] != 'immortalized cell line':
            donor_accession = biosample['donor']['accession']
            biosample_search_url = 'https://www.encodeproject.org/search/?searchTerm=%s&type=biosample&biosample_type=stem cell&format=json&frame=object' % donor_accession
            biosamples_associated_with_donor = get_json_from_encode(biosample_search_url)
            if len(biosamples_associated_with_donor) > 0:
                sample_attribute['lineage'] = biosamples_associated_with_donor[0]['biosample_term_name']
            else:
                print 'No stem cell biosamples associated with donor %s' % donor_accession
        else:
            sample_attribute['lineage'] = ''


def add_SA_primary_cell(sample_attribute, biosample):
    if sample_attribute['biomaterial_type'] == 'Primary Cell':
        sample_attribute['cell_type'] = biosample['biosample_term_name']


def add_SA_primary_tissue(sample_attribute, biosample):
    if sample_attribute['biomaterial_type'] == 'Primary Tissue':
        sample_attribute['tissue_type'] = biosample['biosample_term_id']


def create_experiment_attributes(experiment):
    return {
        "assay_type": experiment['assay_term_name'],
        "experiment_ontology_uri": 'http://purl.obolibrary.org/obo/%s' % experiment['assay_term_id']
    }

############################# Get ENCODE data #############################
# These functions get ENCODE data from the webservices.

def prep_query_url(url, assembly, taxon_id, limit='all'):
    parsed = urlparse.urlparse(url)
    queries = urlparse.parse_qs(parsed.query)
    queries['limit'] = [limit]
    queries['frame'] = ['object']
    queries['format'] = ['json']
    queries['assembly'] = [assembly]
    queries['replicates.library.biosample.donor.organism.taxon_id'] = [str(taxon_id)]
    query = '&'.join(['&'.join([key + '=' + value for value in values]) for key, values in queries.iteritems()])
    new_url = urlparse.urlunparse((parsed.scheme, parsed.netloc, parsed.path, parsed.params, query, parsed.fragment))
    print new_url
    return new_url


def get_json_from_encode(url):
    # Force return from the server in JSON format
    HEADERS = {'accept': 'application/json'}

    # GET the object
    response = requests.get(url, headers=HEADERS)

    # Extract the JSON response as a python dict
    json_response = response.json()

    return json_response if '@graph' not in json_response else json_response['@graph']