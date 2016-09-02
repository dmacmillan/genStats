class Bed:
    
    def __init__(self, chrom, start, stop, name=None, score=None, strand=None, thickStart=None, thickEnd=None, itemRgb=None, blockCount=None, blockSizes=None, blockStarts=None):
        self.chrom = chrom
        self.chromStart = int(start)
        self.chromEnd = int(stop)
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSizes = blockSizes
        self.blockStarts = blockStarts
        if score:
            self.score = int(score)
        self.strand = strand
        if thickStart:
            self.thickStart = int(thickStart)
        if thickEnd:
            self.thickEnd = int(thickEnd)
        self.itemRgb = itemRgb
        if blockCount:
            self.blockCount = int(blockCount)
        if blockSizes:
            self.blockSizes = [int(x) for x in filter(None, blockSizes.split(','))]
        if blockStarts:
            self.blockStarts = [int(x) for x in filter(None, blockStarts.split(','))]

    def __str__(self):
        return '{}\t{}\t{}\t{}\t{}\t{}'.format(self.chrom, self.chromStart, self.chromEnd, self.name, self.score, self.strand)
