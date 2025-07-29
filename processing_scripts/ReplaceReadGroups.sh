#!/usr/bin/env bash

# *******************************************
# Replace SAMPLE and LIBRARY or
# only SAMPLE across read groups in BAM input file.
# Index the ouput file.
# *******************************************

## Usage
USAGE="Usage: ReplaceReadGroups.sh -i <input_file_bam> -s <sample_name> [-l <library>] [-o output_prefix]"

## Default args
output_prefix="output"
nt=$(nproc)

## Functions
check_args()
{
  arg_names=($@)
  for arg_name in ${arg_names[@]}; do
      [ -z ${!arg_name} ] && \
      echo "Missing Argument <${arg_name}>" 1>&2 && \
      echo $USAGE 1>&2 && \
      exit 1
  done
  return 0
}

## Bash command line definition
while getopts 'i:s:l:o:h' opt; do
  case $opt in
    # Required arguments
    i) input_file_bam=${OPTARG} ;;
    s) sample_name=${OPTARG} ;;
    # Optional arguments
    l) library=${OPTARG} ;;
    o) output_prefix=${OPTARG} ;;
    ?|h)
      echo $USAGE 1>&2
      exit 1
      ;;
  esac
done
shift $(($OPTIND -1))

## Check arguments
check_args input_file_bam sample_name

# Update @RG
if [ ! -z ${library} ]; then
  samtools view --no-PG -H $input_file_bam | \
    sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" -e "/^@RG/ s/LB:[^\t]*/LB:${sample_name}.${library}/" | \
    samtools reheader --no-PG - $input_file_bam > ${output_prefix}.bam || exit 1
else
  samtools view --no-PG -H $input_file_bam | \
    sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" | \
    samtools reheader --no-PG - $input_file_bam > ${output_prefix}.bam || exit 1
fi

# Index
samtools index -@ $nt ${output_prefix}.bam || exit 1

# Check BAM integrity
py_script="
import sys, os

def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    fb = open(filename, 'rb')
    fb.seek(size - 28)
    EOF = fb.read(28)
    fb.close()
    if EOF != EOF_hex:
        sys.stderr.write('EOF is missing\n')
        sys.exit(1)
    else:
        sys.stderr.write('EOF is present\n')

check_EOF('${output_prefix}.bam')
"

python -c "$py_script" || exit 1
