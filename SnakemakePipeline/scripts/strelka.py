import os
import re
from Bio import SeqIO

def strelka(refF, healthyF, tumorF, output):
    for record in SeqIO.parse(refF, 'fasta'):
        region = record.id + ':1-' + str(len(record.seq))
        break
    os.system('strelka2 --ref '+ refF + ' --region ' + region + ' --normal-align-file ' + healthyF + ' --tumor-align-file ' + tumorF + ' --somatic-snv-file ' + output)
    print('Executed Strelka successfully')

strelka(str(snakemake.input['ref']), str(snakemake.input['healthy']), str(snakemake.input['tumor']), str(snakemake.output))
