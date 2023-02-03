import pysam


def coverageFromBam(path, allele, variantPos):
    bamFile = pysam.AlignmentFile(path, 'rb')
    cov = 0
    for i in range(4):
        # count coverage returns [4 x len reference] matrix, one row per nucleotide and
        # for each position the number of occurrences per nucleotide
        cov += bamFile.count_coverage(contig=allele)[i][variantPos]
    return cov