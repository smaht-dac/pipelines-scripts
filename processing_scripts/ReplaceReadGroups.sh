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
      echo "Mising Argument <${arg_name}>" 1>&2 && \
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
