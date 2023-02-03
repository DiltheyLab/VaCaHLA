import pysam
import re

def remove_Read_Pairs(input_File, count):
    counter = 0
    keep_reads = {}
    bam = pysam.AlignmentFile(input_File, 'rb')
    for read in bam:
        if read.is_unmapped or read.is_secondary:
            continue
        if read.qname not in keep_reads:
            keep_reads[read.qname] = read
        else:
            ref1 = re.search('HLA-.*\*', read.reference_name).group(0)
            ref2 = re.search('HLA-.*\*', keep_reads[read.qname].reference_name).group(0)
            if ref1 == ref2:
                if read.is_read1:
                    yield read, keep_reads[read.qname]
                else:
                    yield keep_reads[read.qname], read
            else:
                counter += 1
                print('Unequal references! ' + ref1 + ' ' + ref2)
            del keep_reads[read.qname]

    with open(count, 'w') as file:
        file.write(str(counter))


input = snakemake.input['bam']
count = snakemake.output[0]

template = pysam.AlignmentFile(input, 'rb')
output = pysam.AlignmentFile(snakemake.output[1], 'wb', template=template)

for read1, read2 in remove_Read_Pairs(input, count):
    output.write(read1)
    output.write(read2)

output.close()

