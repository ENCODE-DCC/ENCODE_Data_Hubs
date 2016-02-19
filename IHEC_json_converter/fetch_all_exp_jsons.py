import sys
import getopt
import json
from datetime import datetime
import rnaseq, bisulfite, chipseq

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

    #Todo: Merge experiments as a single JSON

    #Whole-Genome Bisulfite Sequencing experiments
    print("Processing WGB-Seq...")
    try:
        data = bisulfite.bisulfite_wrapper(assembly=assembly, taxon_id=taxon_id)
        filename = 'WGB-Seq_%s_%s_%s.json' % (taxon_id, assembly, date_str)
        output_file(data, filename)
        print("Done.")
    except Exception as e:
        print('An error occured while fetching WGB-Seq experiments: ' + e.message)
    print

    #RNA-Sequencing experiments
    print("Processing RNA-Seq...")
    try:
        data = rnaseq.rna_seq_wrapper(assembly=assembly, taxon_id=taxon_id)
        filename = 'RNA-Seq_%s_%s_%s.json' % (taxon_id, assembly, date_str)
        output_file(data, filename)
        print("Done.")
    except Exception as e:
        print('An error occured while fetching RNA-Seq experiments: ' + e.message)
    print

    #ChIP-Seq experiments
    targets = ['H3K27ac', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K4me3', 'H3K9me3']
    for t in targets:
        print("Processing ChIP-Seq %s..." % t)
        try:
            data = chipseq.chip_seq_wrapper(assembly=assembly, taxon_id=taxon_id, target=t)
            filename = 'ChIP-Seq_%s_%s_%s_%s.json' % (taxon_id, assembly, t, date_str)
            output_file(data, filename)
            print("Done.")
        except Exception as e:
            print('An error occured while fetching ChIP-Seq %s experiments: ' % t + e.message)
        print
    print("Operation completed.")



def output_file(data, filename):
    with open('../output/%s' % (filename), 'w+') as outfile:
        json.dump(data, outfile, indent=4)


if __name__ == "__main__":
    main(sys.argv[1:])

