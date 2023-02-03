import os
import re
from Bio import SeqIO

def get_Allele_Names(refFile, bamFile, gene, patient, path, output):
    recordList = []
    with open(refFile, 'r') as fasta:
        records = SeqIO.parse(fasta, 'fasta')
        for record in records:
            if re.search('HLA-' + gene + '\*', record.id):
                recordList.append(record.id)

    if not recordList:
        exit('Warning! No alleles found for gene ' + gene + ' patient ' + patient + ' file ' + refFile)
    else:
        print('samtools view -F 2308 ' + str(bamFile) + ' ' + ' '.join(recordList) + ' -b > ' + str(output))
        os.system('samtools view -F 2308 ' + str(bamFile) + ' ' + ' '.join(recordList) + ' -b > ' + str(output))

    # ensure python script was not interrupted
    os.system('touch ' + path + '/assign_reads_done_' + patient + '_' + gene + '.txt')


get_Allele_Names(snakemake.input['ref'], snakemake.input['bam'], snakemake.wildcards['gene'], str(snakemake.wildcards['patientID']), snakemake.wildcards['path'], snakemake.output[0])
