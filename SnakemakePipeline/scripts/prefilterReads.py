'''
After the first alignment some reads might be mapped outside the HLA region. We only keep those reads that map to some HLA gene.
Here, we filter for those reads + only where both mates are mapped + only primary + no supplementary alignments (Flag -F 2308).
'''

import os
import re
from Bio import SeqIO

def get_allele_names(refFile, bamFile, patient, path, r1_fq, r2_fq):
    recordList = []
    with open(refFile, 'r') as fasta:
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            if re.search('HLA-.+\*', record.id):
                recordList.append(record.id)

    if not recordList:
        exit('Warning! No alleles found for patient ' + patient + ' file ' + refFile)
    else:
        print('samtools view -F 2308 ' + str(bamFile) + ' ' + ' '.join(recordList) + ' -b | samtools sort -n | samtools fastq -1 ' + str(r1_fq) + ' -2 ' + str(r2_fq) + ' -0 /dev/null -s /dev/null -n')
        os.system('samtools view -F 2308 ' + str(bamFile) + ' ' + ' '.join(recordList) + ' -b | samtools sort -n | samtools fastq -1 ' + str(r1_fq) + ' -2 ' + str(r2_fq) + ' -0 /dev/null -s /dev/null -n')

    # ensure python script was not interrupted
    os.system('touch ' + path + '/prefilter_reads_done_' + patient + '.txt')

get_allele_names(snakemake.input['ref'], snakemake.input['bam'], str(snakemake.wildcards['patientID']), snakemake.wildcards['path'], snakemake.output[0], snakemake.output[1])
