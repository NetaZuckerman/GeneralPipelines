#!/bin/bash
###############################
# Dana Bar-Ilan
# 02.08.2021
#
# in general:
# mapping fq files to ref-seq <bwa> -> keep mapped reads, sort, index <samtools> -> create consensus sequence <samtools & bcftools>
# -> align consensuses to ref-seq <mafft>
# + produce report

# requirements: samtools v1.10, bcftools v1.9, ivar v. 1.2.2, mafft v7.215
trap "kill 0" EXIT
eval "$(conda shell.bash hook)"
conda activate CoronaPipeline

path=`dirname "${0}"`

# TODO: keep errors in log file to review later
# TODO: add more documentations

function initialize_globals() {
  dirs_flag=false
  threads=32
  input_path=""
  single_end=false
}

# parse input with flags
function usage() {
    cat <<EOF
Usage: $0 [options]

required: [-h | -d | -i AND -r]
-h| --help                      print this usage message and exit. Ignore the rest
-i              [fastq.gz/path] input path to fastq.gz files location
-r|--refseq     [refseq/path/]  user defined reference. required: refseq/path/ - path to reference fasta file

optional:
--threads       [int]           number of threads. default: 32
-s|--single-end                 single end sequencing

EOF
exit 0
}

function get_user_input() {
  while (( "$#" )); do
    case "$1" in
      -r|--refseq)
        shift
        refseq="$1"
        shift
        ;;
      --threads)
        shift
        threads="$1"
        shift
        ;;
      -i)
        shift
        input_path="$1"
        shift
        ;;
      -h|--help)
        usage
        shift
        ;;
      -s|--single-end)
        single_end=true
        shift
        ;;
      -*|--*=)
        echo "Error: Unsupported flag $1" >&2
        exit 1
        ;;
    esac
  done
}

function check_flags() {
  if [ -z "$input_path" ]; then
    echo "Please provide -i <input_path> to your fastq.gz files location." >&2
    usage
    exit 1
  fi

  if [ -z "$refseq" ]; then
    mkdir -p refs
    echo please provide reference sequence --refseq! >&2
    exit 1
  fi
}


function check_input() {
  # check that input path ends with '/', else add '/' to the end of the string.
    if [[ $input_path != */ ]]; then
      input_path="$input_path"/
    fi
}

function single_end_mapping() {
  bwa index "$refseq"
  if [ -d BAM/ ]; then # not first run, rm BAM files to avoid mixups
    rm BAM/* 2> /dev/null # if file not found it's ok. no need to see on screen
  fi
  if [ -d CNS/ ]; then
    rm CNS/* 2> /dev/null
  fi
  if [ -d CNS_5/ ]; then
    rm CNS_5/* 2> /dev/null
  fi
  mkdir -p BAM CNS alignment results

  for file in "$input_path"*.fastq*; do
    if [[ $file == *Undetermined*.fastq* || $file == *unpaired*.fastq*  || $file == *singletons* ]]; then # ignore undetermined and singletons
        continue
    fi

    if [[ $file == *.gz ]]; then
      output=`basename  "$file" .fastq.gz`
    else
      output=`basename "$file" .fastq`
    fi
    # replace to short name:
    if [[ $output == Sh_* ]]; then
      output=${output/Sh_/}
    fi

    output=BAM/$output.bam

    bwa mem -v1 -t"$threads" "$refseq" "$file" | samtools view -@ "$threads" -b - > $output
  done
}

# paired end mapping: map reads to reference. straight to bam format.
function map_to_ref() {
  # index reference
  bwa index "$refseq"

  if [ -d BAM/ ]; then # not first run, rm existing files to avoid mix ups
    rm BAM/* 2> /dev/null #  if file not found it's ok. no need to print error to screen.
  fi
  if [ -d CNS/ ]; then
    rm CNS/* 2> /dev/null
  fi
  if [ -d CNS_5/ ]; then
    rm CNS_5/* 2> /dev/null
  fi
  mkdir -p BAM CNS alignment results  # -p: only if directory does not exist. else continue without errors.

  for r1 in "$input_path"*R1*.fastq*; do
    if [[ $r1 == *Undetermined*.fastq* || $r1 == *unpaired*.fastq*  || $r1 == *singletons* ]]; then # ignore undetermined and singletons
        continue
    fi
    r2=${r1/R1/R2}
    output=${r1/_R1/}
    output=${output/_paired/}
    if [[ $file == *.gz ]]; then
      output=`basename ${output/.fastq.gz/}`
    else
      output=`basename "$output" .fastq`
    fi
    # replace to short name:
    if [[ $output == Sh_* ]]; then
      output=${output/Sh_/}
    fi
    output=$( echo "$output" | cut -d'_' -f 1 )  # SHORT NAME
    output=BAM/$output.bam

    bwa mem -v1 -t"$threads" "$refseq" "$r1" "$r2" | samtools view -@ "$threads" -b - > $output
  done
}

function keep_mapped_reads() {
  for file in BAM/*.bam; do
     samtools view -@ "$threads" -b -F 260 $file > BAM/`basename $file .bam`.mapped.bam
  done
}

function sort_index_bam() {
  for file in BAM/*.mapped.bam; do
    samtools sort -@ "$threads" "$file" -o BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam
    samtools index -@ "$threads" BAM/`basename "$file" .mapped.bam`.mapped.sorted.bam
  done
}


function depth() {
    mkdir -p QC/depth/
    # create depth files from each mapped sorted and indexed bam file. -a: include all sample's positions including 0 depth.
    for file in BAM/*.mapped.sorted.bam; do
      samtools depth -a "$file" > QC/depth/`basename "$file" .mapped.sorted.bam`.txt
    done
}

function consensus() {
  mkdir -p CNS CNS_5
  for file in BAM/*.mapped.sorted.bam; do
    # ivar instead of bcftools:

    file_name=`basename "$file" .mapped.sorted.bam`
    if [ "$single_end" == FALSE ]; then  # if paired end leave quality threshold 20
      # CNS1
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 1 -p CNS/"$file_name"
      # CNS5
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 5 -p CNS_5/"$file_name"
    else  # of single end allow low quality bases (MinIon)
      # CNS1
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 1 -p CNS/"$file_name" -q 10
      # CNS5
      samtools mpileup -A "$file" | ivar consensus -t 0.6 -m 5 -p CNS_5/"$file_name" -q 10
    fi
  done
  # remove qual files:
  rm CNS/*.qual.txt CNS_5/*.qual.txt
}

# change fasta header from ref to sample name - not needed when using ivar
function change_fasta_header() {
  for file in CNS/*.fa*; do
    # change header to sample name:
    name=`basename $file`
    sed -i "s/>.*/>${name%%.*}/" "$file"
  done

  for file in CNS_5/*.fa*; do
    # change header to sample name:
    name=${file/CNS_5\//} # ${var/find/replace} => remove 'CNS/' prefix
    sed -i "s/>.*/>${name%%.*}/" "$file"
  done
}

function mafft_alignment() {
  # align with MAFFT
  # https://towardsdatascience.com/how-to-perform-sequence-alignment-on-2019-ncov-with-mafft-96c1944da8c6
  conda activate --stack nextstrain

  cat CNS_5/*.fa* > alignment/all_not_aligned.fasta
#  mafft --clustalout alignment/all_not_aligned.fasta > alignment/all_aligned.clustalout
#  mafft alignment/all_not_aligned.fasta > alignment/all_aligned.fasta
  augur align \
  --sequences alignment/all_not_aligned.fasta \
  --reference-sequence "$refseq" \
  --output alignment/all_aligned.fasta

  conda deactivate  # deactivate nextstrain only
}


function results_report() {
  report=QC/report.txt
  # samtools coverage headers: 1#rname  2startpos  3endpos  4numreads  5covbases  6coverage  7meandepth  8meanbaseq  9meanmapq
  echo -e "sample\tmapped%\tmappedreads\ttotreads\tcovbases\tcoverage%\tcoverageCNS_5%\tmeandepth\tmaxdepth\tmindepth" > "$report"
  for file in BAM/*.mapped.sorted.bam; do

    sample_name=`basename $file .mapped.sorted.bam`
    if [[ $sample_name == Sh_* ]]; then
      sample_name=${sample_name/Sh_/}
    fi

    if [ "$single_end" == false ]; then
      sample_name=$( echo "$sample_name" | cut -d'_' -f 1 ) # SHORT NAME
    fi

    original_bam=${file/.mapped.sorted.bam/.bam} # {var/find/replace}
    tot_reads=$(samtools view -c "$original_bam") # do not use -@n when capturing output in variable
    coverage_stats=( $(samtools coverage -H "$file" | cut -f4,5,6) ) # number of mapped reads, covered bases, coverage
    breadth_cns5=$(cut -f3 QC/depth/`basename $file .mapped.sorted.bam`.txt | awk '$1>5{c++} END{print c+0}')
    genome_size=$(cat QC/depth/`basename $file .mapped.sorted.bam`.txt | wc -l)
    coverage_cns5=$(echo "$breadth_cns5/$genome_size"*100 | bc -l)
    mapped_num=${coverage_stats[0]}
    percentage_mapped=$(awk -v m="$mapped_num" -v t="$tot_reads" 'BEGIN {print (m/t)*100}')
    depths=$(awk '{if($3==0){next}; if(min==""){min=max=$3}; if($3>max) {max=$3}; if($3< min) {min=$3}; total+=$3; count+=1} END {print total/count"\t"max"\t"min}' QC/depth/`basename $file .mapped.sorted.bam`.txt)
    echo -e "${sample_name}\t${percentage_mapped}\t${mapped_num}\t${tot_reads}\t${coverage_stats[1]}\t${coverage_stats[2]}\t${coverage_cns5}\t${depths}" >> "$report"
  done
}

########################### MAIN ###############################
# call all functions
# user input:
initialize_globals
get_user_input "$@"
check_flags
check_input
if [ "$single_end" == true ]; then
  single_end_mapping
else
  # start workflow:
  map_to_ref
fi
keep_mapped_reads
sort_index_bam
depth
consensus
if [ "$single_end" == false ]; then
  change_fasta_header
fi

mafft_alignment
results_report

wait

conda deactivate
echo "pipeline finished! (:"

