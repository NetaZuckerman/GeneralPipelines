# GeneralPipelines
General viruses NGS data processing pipelines for alignment, consensus creation, etc.

1. install conda environments:

env/GeneralPipeline.yml
env/nextstrain.yml


2. further dependencies:
samtools V.1.10 and above.
   
3. run:

bash generalPipe.sh -i fastq/dir/path -r refseq/fasta/path

---

### The pipeline:
1. for each fastq pair: map to reference sequence (bwa mem)

2. create consensus sequence for each fastq pair (ivar)

3. align all samples + reference sequence (mafft via augur)

4. create qc report based on coverage and depth (samtools)