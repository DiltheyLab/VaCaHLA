import re
import os

def mutect2(refF, healthyF, tumorF, gene, tid, hid, path, output):
    x = re.search('/(\w+)Samples/hg38BAM/rg_hg38(_\d+_\w+).bam', healthyF)
    rgsmHealthy = x.group(1) + x.group(2)
    y = re.search('/(\w+)Samples/hg38BAM/rg_hg38(_\d+_\w+).bam', tumorF)
    rgsmTumor = y.group(1) + y.group(2)
    os.system('gatk --java-options "-Xmx4g" Mutect2 -R ' + refF + ' -I ' + tumorF + ' -tumor ' + rgsmTumor + ' -I ' + healthyF + ' -normal ' + rgsmHealthy + ' -O ' + output + ' --max-assembly-region-size 500 --max-num-haplotypes-in-population 300 --max-unpruned-variants 200 --max-reads-per-alignment-start 0')
    print('Executed Mutect2 successfully')

    # ensure python script was not interrupted
    os.system('touch ' + path + '/mutect2/mutect2_done_' + tid + '_' + hid + '_' + gene + '.txt')


mutect2(str(snakemake.input['ref']), str(snakemake.input['healthy']), str(snakemake.input['tumor']), snakemake.wildcards['gene'], snakemake.wildcards['patientIDtumor'], snakemake.wildcards['patientIDhealthy'], snakemake.wildcards['path'], str(snakemake.output[0]))

