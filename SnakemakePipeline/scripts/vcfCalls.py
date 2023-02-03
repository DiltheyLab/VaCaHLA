import vcf

# read from the vcf files and return all variant calls as list of records
def getVariants(path):
    records = []
    vcf_reader = vcf.Reader(open(path,'r'))
    for rec in vcf_reader:
        if not rec.FILTER:
            # print(patientID, gene, rec)
            records.append(rec)
    return records