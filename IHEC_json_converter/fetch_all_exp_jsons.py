import sys
import getopt
import json
from datetime import datetime
from reference_epigenome_experiments import collect_experiments

def main(argv):
    opts, args = getopt.getopt(argv, "", ["assembly=", "taxon-id="])

    assembly = ''
    taxon_id = ''
    for opt, arg in opts:
        if opt == '--assembly':
            assembly = arg
        if opt == '--taxon-id':
            taxon_id = arg

    if assembly == '':
        raise Exception("Missing assembly parameter.")
    if taxon_id == '':
        raise Exception("Missing taxon id.")

    date_str = datetime.now().date()

    filename = 'ENCODE.{}.{}.{}.json'.format(taxon_id, assembly, date_str)
    data = collect_experiments()
    output_file(data, filename)



def output_file(data, filename):
    with open('../output/%s' % (filename), 'w+') as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":
    main(sys.argv[1:])

