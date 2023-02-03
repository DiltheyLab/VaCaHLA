import pysam
from collections import defaultdict


# uses alignment positions to create k-mer
def alignmentWithPosArray(read, variantIndexinRead, k):
    # variantIndexinRead = read.get_reference_positions(True).index(callPos)
    kStart = variantIndexinRead - k//2
    kEnd = variantIndexinRead + k//2
    try:
        if read.query_sequence[kStart] and read.query_sequence[kEnd]:
            return read.query_sequence[kStart:kEnd + 1]
    except Exception as e:
        print('ATTENTION! READ IS TOO SHORT ', read)
        print(variantIndexinRead)
        print(read.get_reference_positions(True))
        print(kStart, kEnd)
        exit(e)


# finds the k-mer that carries the variant and occurs most often
def getKmer(path, sample, patientID, gene, rec, k=31):
    allele = rec.CHROM
    callPos = rec.POS   # for deletion its 0 based

    if len(rec.ALT) > 1:
        print('ATTENTION! MULTIPLE ALTERNATIVES:', rec.ALT, sample, patientID, gene)
        exit()

    ALT = rec.ALT[0]

    callPos -= 1    # vcf starts counting at 1 but has to start at 0
    samfile = pysam.AlignmentFile(path + '/' + sample + '/hg38BAM/rg_hg38_' + patientID + '_' + gene + '.bam', 'rb')

    if sample == 'healthySamples':
        samfile = pysam.AlignmentFile(path + '/' + sample + '/hg38BAM/hg38_bwa_' + patientID + '_' + gene + '.bam', 'rb')

    # stores all variant carrying k-mers together with their number of occurrences
    kmers = defaultdict(lambda: 0)
    # stores variant carrying k-mers together with number of high-quality reads from which k-mer originates
    readQual = defaultdict(lambda: 0)
    # stores variant carrying k-mers together with base quality of variant
    baseQual = defaultdict(lambda: 0)
    # variantReadsCount = [number of variant carrying reads, number of + reads, number of - reads]
    variantReadsCount = [0, 0, 0]

    # get all reads covering the variant position
    # .fetch last position is exclusive, additionally it uses zero based values
    for read in samfile.fetch(allele, callPos, callPos + 1):
        if read.is_secondary or read.is_supplementary or read.is_unmapped:
            continue
        
        # get affected_end of variant, notice that it is exclusive for deletions and snps (?)
        variantEnd = rec.affected_end

        # check whether read covers variant position
        try:
            variantIndexinRead = read.get_reference_positions(True).index(callPos)
            indexVariantEnd = variantIndexinRead + 1 if rec.is_snp else read.get_reference_positions(True).index(variantEnd)
        except Exception as e:
            # print(e)
            continue

        if rec.is_snp and len(read.query_sequence[variantIndexinRead:indexVariantEnd]) > 1:
            print('is snp but variant is too long!!! check it out!!!')
            exit()

        # if read does not carry variant continue
        if read.query_sequence[variantIndexinRead:indexVariantEnd] != ALT:
            continue

        variantReadsCount[0] += 1
        if read.is_reverse:
            variantReadsCount[2] += 1
        else:
            variantReadsCount[1] += 1

        kStart = variantIndexinRead - ((k + 1) // 2) - len(ALT)
        kEnd = variantIndexinRead + ((k + 1) // 2) + len(ALT)
        readLength = len(read.get_reference_positions(True))

        # get only those reads that cover region completely
        if kStart > 0 and kEnd < readLength:

            if rec.is_snp:
                kmer = alignmentWithPosArray(read, variantIndexinRead, k)
                if not len(kmer) == 31:
                    print('ATTENTION! SNP kmer is not 31 bases long! Length: ', len(kmer), kmer)
                    print('Patient; Gene; Allele; Call position; Variant: ', patientID, gene, allele, callPos, ALT)
                    exit(1)

                kmers[kmer] += 1  # count occurrences of each k-mer
                # count the number of high-quality variant alleles
                if read.mapq >= 30:
                    readQual[kmer] += 1
                # count the number of hig base qualities
                if read.query_qualities[variantIndexinRead] >= 30:
                    baseQual[kmer] += 1

            elif rec.is_indel:
                # call adds one normal mapped base, we keep it here as it is part of the k-mer in any case
                indelLength = len(ALT)

                # kmerCenterPosInRead is the start position from which we will construct the k-mer k//2 bases to the left and k//2 to the right
                kmerCenterPosInRead = variantIndexinRead + indelLength // 2
                kmer = alignmentWithPosArray(read, kmerCenterPosInRead, k)
                if not len(kmer) == k:
                    print('ATTENTION! Indel kmer is not 31 bases long! Length: ', len(kmer), kmer)
                    print('Patient; Gene; Allele; Call position; Variant: ', patientID, gene, allele, callPos, ALT)
                    exit(1)

                kmers[kmer] += 1  # count occurrences of each k-mer
                # count the number of high-quality variant alleles
                if read.mapq >= 30:
                    readQual[kmer] += 1
                # count the number of hig base qualities, deletions automatically have a high base quality, as there is no base to check
                if rec.is_deletion:
                    baseQual[kmer] += 1
                elif read.query_qualities[variantIndexinRead] >= 30:
                    baseQual[kmer] += 1

            else:
                print('What variant type is this?', rec.var_type)
                exit()

        else:
            continue

    samfile.close()
    if not kmers:
        #sys.exit('ATTENTION! List of KMERs is NONE!')
        return None, 0, variantReadsCount, 0, 0

    #totalVariantCoverage = sum(kmers.values())
    maxCountedKmer = max(kmers, key=kmers.get)
    return maxCountedKmer, kmers[maxCountedKmer], variantReadsCount, readQual[maxCountedKmer], baseQual[maxCountedKmer]
