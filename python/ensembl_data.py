__author__ = 'Patrick Short'

# dependencies #
from bioservices.ensembl import Ensembl as e

# Load well-covered regulatory regions from DDD 4k exome data
unqueryable = []
with open("../data/DDD_TRRs.annotated.highcov.txt", 'r') as wc_regs:
    with open("../data/DDD_TRRs.annotated.highcov.sequence.txt", 'w') as out_file:
        header = wc_regs.readline().rstrip().split("\t")
        header.append("seq")
        new_header = "\t".join(header)
        print>>out_file, new_header
        for line in wc_regs:
            line = line.rstrip().split("\t")  # first col is chr, second is start, third is stop
            chr = line[0]
            start = line[1]
            stop = line[2]
            query = chr + ":" + start + ".." + stop + ":1"  # e.g. X:1000000..1000100:1

            try:
                seq = e().get_sequence_by_region(query, species="human", coord_system_version="GRCh37.p13").get('seq')
            except AttributeError:  # seq is integer 400 - means server error. this occurred at end of chr 4 (off edge)
                print 'Could not query region...'
                unqueryable.append(line)
                continue

            if seq:
                line.append(seq.upper())
                print>>out_file, "\t".join(line)
            else:
                print "Could not query region..."
                unqueryable.append(line)
