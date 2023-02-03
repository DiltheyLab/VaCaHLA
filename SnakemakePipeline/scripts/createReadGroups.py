import os
import re

def create_Read_Groups(inputF, gene, patient, path, output):
    x = re.search('/(\w+)Samples/hg38BAM/rg_hg38(_\d+_\w+).bam', str(output))
    if not x:
    	print('Error! Could not match /(\w+)Samples/hg38BAM/rg_hg38(_\d+_\w+).bam')
    	exit(-1)
    	
    rgsm = x.group(1) + x.group(2)

    os.system('picard AddOrReplaceReadGroups I=' + str(inputF) + ' O=' + str(output) + ' RGID=4  RGLB=lib1 RGPL=illumina  RGPU=unit1  RGSM=' + rgsm)
    os.system('samtools index ' + str(output))
    #print('Executed createReadGroups successfully')

    # ensure python script was not interrupted
    os.system('touch ' + path + '/hg38BAM/rg_done_' + patient + '_' + gene + '.txt')


create_Read_Groups(snakemake.input, snakemake.wildcards['gene'], str(snakemake.wildcards['patientID']), snakemake.wildcards['path'], snakemake.output[0])
