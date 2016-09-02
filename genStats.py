import argparse, os, sys, math
from Bed import *
import pysam
import numpy as np

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
parser.add_argument('-g', '--gene', nargs= '+', default='result', help='Only look at a specific set of genes')
parser.add_argument('-n', '--name', default='result', help='Name for the file output. Default is "result"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Directory to output to. Will be created if it does not exist. Default is "{}"'.format(os.getcwd()))

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

outfile = os.path.join(args.outdir, args.name + '.stats')

def genStats(alignment, regions, out=None, name='NA', cutoff=100):
    aln = pysam.AlignmentFile(alignment, 'rc')
    header = ('\t').join(['CHROM', 'GENE', 'STRAND',
                          'SAMPLE', 'START', 'END',
                          'LENGTH', 'TOTAL', 'AVERAGE',
                          'MIN', 'Q1', 'MEDIAN', 'Q3', 
                          'MAX', 'SEM'])
    with open(out, 'w') as o:
        o.write(header + '\n')
        with open(regions, 'r') as f:
            for line in f:
                bed = Bed(*line.strip().split('\t'))
                pileup = []
                for nuc in aln.pileup(bed.chrom, bed.chromStart, bed.chromEnd, truncate=True):
                    #if bed.chromStart <= nuc.pos <= bed.chromEnd:
                    pileup.append(nuc.n)
                if not pileup:
                    continue
                _length = len(pileup)
                _total = sum(pileup)
                _average = _total/_length
                if _average < cutoff:
                    continue
                _median = np.percentile(pileup, 50)#median(pileup)
                _q1 = np.percentile(pileup, 25)
                _q3 = np.percentile(pileup, 75)
                _max = max(pileup)
                _min = min(pileup)
                _sem = np.std(pileup)/math.sqrt(_length)
                line = ('\t').join([str(x) for x in [bed.chrom, 
                                                     bed.name,
                                                     bed.strand,
                                                     name,
                                                     bed.chromStart,
                                                     bed.chromEnd,
                                                     _length,
                                                     _total,
                                                     _average,
                                                     _min,
                                                     _q1,
                                                     _median,
                                                     _q3,
                                                     _max,
                                                     _sem]])
                o.write(line + '\n')

print 'Writing statistics to: {} ...'.format(outfile)
genStats(args.alignment, args.regions, name=args.name, out=outfile)
