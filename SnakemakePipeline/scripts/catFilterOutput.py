# create key from specific columns that have to match between the same variant calls
def create_key(list, file):
    try:
        return (list[0], list[1], list[2], list[4], list[6], list[7])
    except Exception as e:
        print(e, line, file)


mutect1 = snakemake.input['mutect1']
mutect2 = snakemake.input['mutect2']
strelka2 = snakemake.input['strelka2']
out = snakemake.output[0]

print(out)

filterMap = {}

with open(mutect1, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        l = line.strip().split(',')
        key = create_key(l, mutect1)
        filterMap[key] = l + ['1', '0', '0']

with open(mutect2, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        l = line.strip().split(',')
        key = create_key(l, mutect2)
        if key in filterMap:
            filterMap[key][-2] = '1'
        else:
            filterMap[key] = l + ['0', '1', '0']

with open(strelka2, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        l = line.strip().split(',')
        key = create_key(l, strelka2)
        if key in filterMap:
            filterMap[key][-1] = '1'
        else:
            filterMap[key] = l + ['0', '0', '1']

with open(out, 'w') as f:
    f.write('#Kmer,ID Normal,ID Tumor,Is Filtered,Gene,Call Type,Variant,Pos,T alignment coverage,T jf count,T kmer count,T total variant reads,Variant Frequency,Strand bias,T + variant reads,T - variant reads,T high-quality reads,T high base quality,N coverage,N jf count,N variant reads,N + variant reads,N - variant reads\n')
    for line in filterMap.values():
        f.write(','.join(line) + '\n')
