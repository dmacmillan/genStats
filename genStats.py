import argparse, os, sys, math
from Bed import *
import pysam
import numpy as np
from multiprocessing import Pool
from functools import partial

def median(lst):
    sortedLst = sorted(lst)
    lstLen = len(lst)
    index = (lstLen - 1) // 2

    if (lstLen % 2):
        return sortedLst[index]
    else:
        return (sortedLst[index] + sortedLst[index + 1])/2.0

parser = argparse.ArgumentParser(description='')

parser.add_argument('alignment', help='Alignment file in CRAM format')
parser.add_argument('regions', help='Bed file of regions to generate statistics for')
parser.add_argument('-p', '--processes', type=int, default=1, help='Number of processes to use')
parser.add_argument('-g', '--gene', nargs= '+', default='result', help='Only look at a specific set of genes')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output to. Will be created if it does not exist. Default is "{}"'.format(os.getcwd()))

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

outfile = os.path.join(args.outdir, args.name + '.stats')

# alignment_file: a string representing the path to an alignment file in cram format
# region: a Bed object
def genStatForRegion(region, alignment_file, name='NA'):
    aln = pysam.AlignmentFile(alignment_file, 'rc')
    pileup = [nuc.nsegments for nuc in aln.pileup(region.chrom, region.chromStart, region.chromEnd, truncate=True)]
    if not pileup:
        return None
    _length = len(pileup)
    _total = sum(pileup)
    _average = _total/_length
    _median = np.percentile(pileup, 50)
    _q1 = np.percentile(pileup, 25)
    _q3 = np.percentile(pileup, 75)
    _max = max(pileup)
    _min = min(pileup)
    _sem = np.std(pileup)/math.sqrt(_length)
    return ('\t').join([str(x) for x in [region.chrom, 
                                         region.name,
                                         region.strand,
                                         name,
                                         region.chromStart,
                                         region.chromEnd,
                                         _length,
                                         _total,
                                         _average,
                                         _min,
                                         _q1,
                                         _median,
                                         _q3,
                                         _max,
                                         _sem]])

def parseBed(bedfile):
    results = []
    with open(bedfile, 'r') as f:
        for line in f:
            bed = Bed(*line.strip().split('\t'))
            results.append(bed)
    return results

if __name__ == '__main__':
    header = ('\t').join(['CHROM', 'GENE', 'STRAND',
                          'SAMPLE', 'START', 'END',
                          'LENGTH', 'TOTAL', 'AVERAGE',
                          'MIN', 'Q1', 'MEDIAN', 'Q3', 
                          'MAX', 'SEM'])
    regions = parseBed(args.regions)
    output = open(outfile, 'w')
    pool = Pool(processes = args.processes)

    partial_genStatForRegion = partial(genStatForRegion, alignment_file=args.alignment, name=args.name)

    lines = pool.map(partial_genStatForRegion, regions)

#def genStats(alignment, regions, out=None, name='NA', cutoff=100):
#    aln = pysam.AlignmentFile(alignment, 'rc')
#    header = ('\t').join(['CHROM', 'GENE', 'STRAND',
#                          'SAMPLE', 'START', 'END',
#                          'LENGTH', 'TOTAL', 'AVERAGE',
#                          'MIN', 'Q1', 'MEDIAN', 'Q3', 
#                          'MAX', 'SEM'])
#    with open(out, 'w') as o:
#        o.write(header + '\n')
#        with open(regions, 'r') as f:
#            for line in f:
#                bed = Bed(*line.strip().split('\t'))
#                pileup = []
#                for nuc in aln.pileup(bed.chrom, bed.chromStart, bed.chromEnd, truncate=True):
#                    #if bed.chromStart <= nuc.pos <= bed.chromEnd:
#                    pileup.append(nuc.n)
#                if not pileup:
#                    continue
#                _length = len(pileup)
#                _total = sum(pileup)
#                _average = _total/_length
#                if _average < cutoff:
#                    continue
#                _median = np.percentile(pileup, 50)#median(pileup)
#                _q1 = np.percentile(pileup, 25)
#                _q3 = np.percentile(pileup, 75)
#                _max = max(pileup)
#                _min = min(pileup)
#                _sem = np.std(pileup)/math.sqrt(_length)
#                line = ('\t').join([str(x) for x in [bed.chrom, 
#                                                     bed.name,
#                                                     bed.strand,
#                                                     name,
#                                                     bed.chromStart,
#                                                     bed.chromEnd,
#                                                     _length,
#                                                     _total,
#                                                     _average,
#                                                     _min,
#                                                     _q1,
#                                                     _median,
#                                                     _q3,
#                                                     _max,
#                                                     _sem]])
#                o.write(line + '\n')
#
#print 'Writing statistics to: {} ...'.format(outfile)
#genStats(args.alignment, args.regions, name=args.name, out=outfile)
