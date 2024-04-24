import sys
import kmerFromReads
import vcfCalls
import variantPosCoverage
import re
from Bio import SeqIO


# create reverse complement if it is lexicographically smaller than k-mer
def reversedComplement(kmer):
    tab = str.maketrans('ACTG', 'TGAC')
    if not kmer < kmer.translate(tab)[::-1]:
        kmer = kmer.translate(tab)[::-1]
    return kmer


# kmer lookup in jellyfish files, returns count of kmer in original reads
def getJellyfishCoverage(kmer, path, pIDtumor, pIDhealthy, gene):
    kmerTumorCount = 0
    kmerNormalCount = 0
    
    if 'N' in kmer:
        print('Kmer Contains N, automatically reject call\nKmer: ' + kmer + ' patient t: ' + pIDtumor + ' h: ' + pIDhealthy + ' gene: ' + gene)
        return kmerTumorCount, 32202

    kmer = reversedComplement(kmer)

    with open(path + '/healthySamples/' + pIDhealthy + '_counts_dumps.fasta','r') as jellyfishHandle:
        for rec in SeqIO.parse(jellyfishHandle, 'fasta'):
            if rec.seq == kmer:
                kmerNormalCount = int(rec.id)
                #print('Healthy sample; Record: ', rec.seq, 'Kmer: ', kmer, rec.id)
                break

    with open(path + '/tumorSamples/' + pIDtumor + '_counts_dumps.fasta','r') as jellyfishHandle:
        for rec in SeqIO.parse(jellyfishHandle, 'fasta'):
            if rec.seq == kmer:
                kmerTumorCount = int(rec.id)
                #print('Tumor sample; Record: ', rec.seq, 'Kmer: ', kmer, rec.id)
                break

    if kmerNormalCount == 0 and kmerTumorCount == 0:
        exit('ATTENTION! Kmer ' + kmer + ' does not exist in jellyfish counts of patient t: ' + pIDtumor + ' h: ' + pIDhealthy + ' gene: ' + gene)

    return kmerTumorCount, kmerNormalCount


# START
path = snakemake.wildcards.path
vcfInput = snakemake.input['vcfInput']
gene = snakemake.wildcards.gene
t = snakemake.wildcards.patientIDtumor
h = snakemake.wildcards.patientIDhealthy
output = str(snakemake.output)


filteredVariants = open(output, 'w')
filteredVariants.write('#Kmer,ID Normal,ID Tumor,Is Filtered,Gene,Call Type,Variant,Pos,T alignment coverage,T jf count,T kmer count,T total variant reads,Variant Frequency,Strand bias,T + variant reads,T - variant reads,T high-quality reads,T high base quality,N coverage,N jf count,N variant reads,N + variant reads,N - variant reads\n')

records = vcfCalls.getVariants(str(vcfInput))
# if no variants called, continue
if not records:
    exit()

# get 31-mer from reads with mutation at pos 15
for rec in records:
    print(rec)
    # kmer filtering snp variants
    if rec.is_snp:
        # what to do if multiple alternative bases for SNP?
        if len(rec.ALT) > 1:
            exit('Multiple alternate variants for SNP: ', rec.ALT)

        tumorKmer, tCount, variantCarryingReadCountsT, qualityCountT, baseQualityCountT = kmerFromReads.getKmer(path, 'tumorSamples', t, gene, rec)
        normalKmer, nCount, variantCarryingReadCountsN, qualityCountN, baseQualityCountN = kmerFromReads.getKmer(path, 'healthySamples', h, gene, rec)

        jellyfishKmerCountTumor = 0
        jellyfishKmerCountNormal = 0
        if tumorKmer:
            jellyfishKmerCountTumor, jellyfishKmerCountNormal = getJellyfishCoverage(tumorKmer, path, t, h, gene)
        else:
            continue

        bamCovAtVariantPos = variantPosCoverage.coverageFromBam(path + '/tumorSamples/hg38BAM/rg_hg38_' + t + '_' + gene + '.bam', rec.CHROM, rec.POS)
        bamCovAtNormalPos = variantPosCoverage.coverageFromBam(path + '/healthySamples/hg38BAM/rg_hg38_' + h + '_' + gene + '.bam', rec.CHROM, rec.POS)

    # MNPs also count as indel here
    elif rec.is_indel:
        # what to do if multiple alternative bases for INDEL?
        if len(rec.ALT) > 1:
            print('Multiple alternate variants for INDEL: ', rec.ALT)
            exit(-1)

        tumorKmer, tCount, variantCarryingReadCountsT, qualityCountT, baseQualityCountT = kmerFromReads.getKmer(path, 'tumorSamples', t, gene, rec)
        normalKmer, nCount, variantCarryingReadCountsN, qualityCountN, baseQualityCountN = kmerFromReads.getKmer(path, 'healthySamples', h, gene, rec)

        jellyfishKmerCountTumor = 0
        jellyfishKmerCountNormal = 0
        if tumorKmer:
            jellyfishKmerCountTumor, jellyfishKmerCountNormal = getJellyfishCoverage(tumorKmer, path, t, h, gene)
        else:
            continue

        bamCovAtVariantPos = variantPosCoverage.coverageFromBam(path + '/tumorSamples/hg38BAM/rg_hg38_' + t + '_' + gene + '.bam', rec.CHROM, rec.POS)
        bamCovAtNormalPos = variantPosCoverage.coverageFromBam(path + '/healthySamples/hg38BAM/rg_hg38_' + h + '_' + gene + '.bam', rec.CHROM,rec.POS)

    else:
        print('Type is unknown, what to do?', t, h, gene, rec)
        exit(-1)

    # add percent values to relative read count values
    percentValuesT = []
    percentValuesN = []
    biasFilterT = []
    for i in range(1, 3):
        p = 0
        if variantCarryingReadCountsT[0] != 0:
            p = (variantCarryingReadCountsT[i]/variantCarryingReadCountsT[0]) * 100
        biasFilterT.append(p)
        percentValuesT.append(str(variantCarryingReadCountsT[i]) + '/' + '{:.0f}'.format(p) + '%')

        p = 0
        if variantCarryingReadCountsN[0] != 0:
            p = (variantCarryingReadCountsN[i]/variantCarryingReadCountsN[0]) * 100
        percentValuesN.append(str(variantCarryingReadCountsN[i]) + '/' + '{:.0f}'.format(p) + '%')

    # add variant frequency
    variantFrequency = round(variantCarryingReadCountsT[0]/bamCovAtVariantPos, 2)

    # Filter variant calls based on kmer counts
    isFiltered = False
    # Global Kmer Filter: if variant kmer count in normal >= 3: reject call (based on jellyfish files)
    if jellyfishKmerCountNormal >= 3:
        isFiltered = True

    # Pileup Frequency Filter
    if variantCarryingReadCountsN[0] >= 2:
        isFiltered = True

    # Quantity and Quality Filter: (currently only for real data)
    if not re.search('testData', path):
        if qualityCountT <= 3:
            isFiltered = True
        if baseQualityCountT <= 3:
            isFiltered = True

    ############################################################################################
    # hard filters
    # if normal coverage too low to compare for variants, reject call
    if  bamCovAtNormalPos < 30:
        isFiltered = True

    # if realtion between variant carrying reads in normal and tumor is high, reject call
    relation = (variantCarryingReadCountsN[0] / variantCarryingReadCountsT[0]) * 100
    if relation > 10:
        isFiltered = True

    strandBias = 'No'
    # filter for strand bias at 25% to 75% relation
    if abs(biasFilterT[0]-biasFilterT[1]) >= 50:
        strandBias = 'Yes'

    ############################################################################################

    filteredVariants.write(
        '{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(
            tumorKmer,
            h,
            t,
            isFiltered,
            gene,
            rec.var_subtype,
            rec.ALT,
            rec.POS,
            bamCovAtVariantPos,
            jellyfishKmerCountTumor,
            tCount,
            variantCarryingReadCountsT[0],
            variantFrequency,
            strandBias,
            percentValuesT[0],
            percentValuesT[1],
            qualityCountT,
            baseQualityCountT,
            bamCovAtNormalPos,
            jellyfishKmerCountNormal,
            variantCarryingReadCountsN[0],
            percentValuesN[0],
            percentValuesN[1]
        )
    )

filteredVariants.close()
