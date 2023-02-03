import re
import sys


read_ids = {}
for l, line in enumerate(open(snakemake.input['r1'])):
    if l % 4 == 0:
        line = line.strip()
        line = re.sub(r' .*', '', line)
        if not line in read_ids:
            read_ids[line] = 0
        read_ids[line] += 1
        if read_ids[line] > 1:
            sys.exit('Read ID occurred multiple times path: ' + file)
