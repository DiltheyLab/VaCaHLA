# ![Logo](vacahla_logo.png)aCaHLA - Variant Calling in the HLA region on Illumina Amplicon Sequencing Data


![Pipeline](SnakemakePipeline.png)

VaCaHLA is designed to perform variant calling int the HLA region using Illumina amplicon sequencing data. The workflow is shown above.
## Pipeline Execution
### Data
As input, the pipeline requires paired-end Illumina amplicon fastq files, personalized reference fasta files containing the individual HLA references, and the GRCh38 reference. A generated test sample and necessary reference data can be found [here]().

### Preparation
This pipeline is based on the Snakemake workflow and management system and uses conda to automatically deploy the required software. Make sure to have [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) installed. Start with creating an environment for the pipeline execution by running the following commands:

```
conda env create -f environment.yaml
conda activate vacahla
```
We use three variant callers, amongst others MuTect, which we run in a singularity container. This must be considered when running the pipeline. Make sure to have singularity installed (requires Linux or Linux VM).

Use the [config file](SnakemakePipeline/envs/config.yaml) to adjust parameters for your execution. The pipeline requires normal and tumor fastq samples and reference sequences as input.  It follows an explanation of the input parameters:

- **prefix**: used as prefix of personal reference fasta files
- **idMapping**: tab separated file listing the tumor ID and normal ID of the read files for each patient in a new line
- **path**: path for input and output;  must contain the folders "healthySamples" and "tumorSamples", containing the normal and tumor fastq files and a subfolder "personalReferences" for the normal data and tumor data, respectively
- **referenceFolder**: subfolder in the previous folders containing personalized references
- **hg38path**: path to the reference sequences contained in GRCh38 of the genes under investigation
- **hg38replaced**: Modified GRCh38 sequence where all target genes (as well as an additional 2000 bases of padding to either side of each gene) is masked with ‘N’ characters and for each target gene, the sequences from the IMGT/HLA database are included
- **sampleGenes**: Genes under investigation
- **prefiltering**: After mapping the reads to personalized references, we filter them. This parameter is used to direct the execution flow of snakemake, such that after filtering the reads are checked for correctness once again.

### Execution
To execute the pipeline use the following command: 

```
snakemake --cores $number$ --use-conda --use-singularity --singularity-args "--bind $HOME"
```

All software requirements are deployed automatically via conda environments.

As output, you get the filtered variant calls as csv files. *False_merged.csv* contains all variant calls from all three variant calling tools that remained after applying our kmer filter. *True_merged.csv* contains the variant calls that were filtered out.
