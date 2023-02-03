import pysam
import re
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from collections import Counter


def count_clipping(file):
    all_counter = 0
    total_clipped = 0
    clip_counts = []
    samfile = pysam.AlignmentFile(file, 'rb')
    for read in samfile.fetch():
        all_counter += 1
        if read.is_unmapped or read.is_supplementary:
            continue
        cigar = read.cigarstring
        if 'S' in cigar:
            total_clipped += 1
            clipping_length = re.findall(r'(\d+)S', cigar)
            if clipping_length:
                clip_counts.extend([int(val) for val in clipping_length])
        # else:
        #    clip_counts.append(0)

    return all_counter, total_clipped, clip_counts


# START
genes = ['A', 'B', 'B2M', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1', 'E', 'F', 'G', 'H', 'TAP1', 'TAP2', 'MICA', 'MICB']
patient = snakemake.wildcards.patientID
path = snakemake.wildcards.path
output = str(snakemake.output)

# Initialization
count_reads = {}
count_clipped = {}
clipping_values = {}

for gene in genes:
    count_reads[gene] = 0
    count_clipped[gene] = 0
    clipping_values[gene] = []

# Get read count information for each gene
for gene in genes:
    file = path + '/hg38BAM/hg38_bwa_' + patient + '_' + gene + '.bam'

    c_reads, c_clipped, clipping = count_clipping(file)

    count_reads[gene] += c_reads
    count_clipped[gene] += c_clipped
    clipping_values[gene].extend(clipping)

# Y labels are used for legend position in plots
y_label_coordinates = {}
for gene in genes:
    val = Counter(clipping_values[gene]).most_common(1)
    #print(val)
    if not val:
        y_label_coordinates[gene] = 0
        continue
    y_label_coordinates[gene] = val[0][1]

# Create histogram plots
fig = make_subplots(rows=8, cols=2, subplot_titles=tuple(genes), vertical_spacing=0.04)

for i in range(0, 16, 2):
    fig.add_trace(go.Histogram(x=clipping_values[genes[i]],
                               name=genes[i],
                               xbins_size=1),
                  row=i // 2 + 1, col=1)
    fig.update_xaxes(title_text='Soft Clipping Length', row=i // 2 + 1, col=1)
    fig.update_yaxes(title_text='Frequency', row=i // 2 + 1, col=1)
    fig.add_annotation(
        text='Total number of reads: ' + str(count_reads[genes[i]]) + '<br>Total number of clipped reads: ' + str(count_clipped[genes[i]]),
        row=i // 2 + 1, col=1, x=200, y=y_label_coordinates[genes[i]], showarrow=False)

    fig.add_trace(go.Histogram(x=clipping_values[genes[i + 1]],
                               name=genes[i + 1],
                               xbins_size=1),
                  row=i // 2 + 1, col=2)
    fig.update_xaxes(title_text='Soft Clipping Length', row=i // 2 + 1, col=2)
    fig.update_yaxes(title_text='Frequency', row=i // 2 + 1, col=2)
    fig.add_annotation(
        text='Total number of reads: ' + str(count_reads[genes[i + 1]]) + '<br>Total number of clipped reads: ' + str(count_clipped[genes[i+1]]),
        row=i // 2 + 1, col=2, x=200, y=y_label_coordinates[genes[i + 1]], showarrow=False)

fig.update_traces(marker=dict(color='midnightblue'), selector=dict(type='histogram'))
fig.update_layout(height=3000, width=1500)
#fig.show()

# Save as html to keep hover information
fig.write_html(output)