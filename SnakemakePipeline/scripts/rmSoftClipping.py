import pysam
import re
import os


def get_clipping(read_list):
    clipping = False
    for read in read_list:
        if not read.cigarstring:
            continue
        if 'S' in read.cigarstring:
            clipping = True
            # clippingLength = max([int(i) for i in re.findall(r'(\d+)S', read.cigarstring)])
            # if clippingLength > 10:
            #     clipping = True
    return clipping


def read_bam(file):
    reads = {}
    is_first, is_sec = False, False
    bamfile = pysam.AlignmentFile(file, 'rb')
    for read in bamfile.fetch():
        if read.is_unmapped or read.is_supplementary:
            continue
        if read.query_name not in reads:
            is_first = read.is_read1
            is_sec = read.is_read2
            reads[read.query_name] = [read]
        else:
            # make sure to get only read pairs! So no two times read 1 or read 2
            if read.is_read1 and is_first:
                continue
            elif read.is_read2 and is_sec:
                continue

            reads[read.query_name].append(read)
            # check whether one of the read mates is soft clipped; if so: remove both, else: return both
            if not get_clipping(reads[read.query_name]):
                yield reads[read.query_name]
            del reads[read.query_name]
        is_first, is_sec = False, False


file = str(snakemake.input['bamfile'])
print(file)

template = pysam.AlignmentFile(file, 'rb')
tmp = str(snakemake.output) + '.tmp'
output = pysam.AlignmentFile(tmp, 'wb', template=template)
print(tmp)
#output = pysam.AlignmentFile('test.bam', 'wb', template=template)

for read_x, read_y in read_bam(file):
    output.write(read_x)
    output.write(read_y)

output.close()

# sort bam file for further processes
os.system('samtools sort ' + tmp + ' -o ' + str(snakemake.output))
os.system('rm ' + tmp)
