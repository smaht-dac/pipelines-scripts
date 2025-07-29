#!/usr/bin/env bash

# *******************************************
# Replace sample name (and optionally library) in all @RG lines
# of a BAM or CRAM file. Reheader and index the output.
# *******************************************

## Usage
usage() {
  cat << EOF
Usage: ReplaceReadGroups.sh [options]

Required arguments:
  -i | --input <input_file>        Input BAM or CRAM file
  -s | --sample <sample_name>      Sample name to replace in all @RG entries

Optional arguments:
  -l | --library <library>         Library name to replace in all @RG entries
  -o | --output <output_prefix>    Output prefix (default: 'output')
  -r | --reference <fasta>         Reference FASTA (required for CRAM input)
  --in_place                       Overwrite original CRAM file in place (CRAM only)
  -h | --help                      Show this help message
EOF
  exit 1
}

## Default args
output_prefix="output"
nt=$(nproc)
in_place=false
is_cram=false

## Functions
err() { echo "[ERROR] $*" >&2; }

check_args() {
  arg_names=("$@")
  for arg_name in "${arg_names[@]}"; do
    if [ -z "${!arg_name}" ]; then
      err "Missing argument: <$arg_name>"
      usage
    fi
  done
}

## Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)
      input_file="$2"
      shift 2
      ;;
    -s|--sample)
      sample_name="$2"
      shift 2
      ;;
    -l|--library)
      library="$2"
      shift 2
      ;;
    -o|--output)
      output_prefix="$2"
      shift 2
      ;;
    -r|--reference)
      reference_fasta="$2"
      shift 2
      ;;
    --in_place)
      in_place=true
      shift
      ;;
    -h|--help)
      usage
      ;;
    *)
      err "Unknown argument: $1"
      usage
      ;;
  esac
done

## Check required args
check_args input_file sample_name

## Detect file extension
input_ext=$(basename "$input_file" | awk -F. '{print tolower($NF)}')
[[ "$input_ext" == "cram" ]] && is_cram=true

## Check CRAM-specific args
if $is_cram; then
  check_args reference_fasta
  if $in_place; then
    echo "[WARNING] Overwriting original CRAM file: $input_file"
  fi
elif $in_place; then
  err "--in_place is only supported for CRAM input files"
  usage
fi

## Modify read group headers
if ! $is_cram; then
  # BAM case
  if [ -n "$library" ]; then
    samtools view --no-PG -H "$input_file" | \
      sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" \
          -e "/^@RG/ s/LB:[^\t]*/LB:${sample_name}.${library}/" | \
      samtools reheader --no-PG - "$input_file" > "${output_prefix}.bam" || exit 1
  else
    samtools view --no-PG -H "$input_file" | \
      sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" | \
      samtools reheader --no-PG - "$input_file" > "${output_prefix}.bam" || exit 1
  fi
else
  # CRAM case
  trap 'rm -f tmp_header' EXIT

  # Create temporary header with modified read groups
  if [ -n "$library" ]; then
    samtools view --no-PG -H -T "$reference_fasta" "$input_file" | \
      sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" \
          -e "/^@RG/ s/LB:[^\t]*/LB:${sample_name}.${library}/" > tmp_header || exit 1
  else
    samtools view --no-PG -H -T "$reference_fasta" "$input_file" | \
      sed -e "/^@RG/ s/SM:[^\t]*/SM:${sample_name}/" > tmp_header || exit 1
  fi

  # Apply header to CRAM
  samtools reheader --no-PG -i tmp_header "$input_file" || exit 1

  # If not in-place, convert to BAM
  if ! $in_place; then
    samtools view -@ "$nt" --no-PG -h --bam -T "$reference_fasta" -o "${output_prefix}.bam" "$input_file" || exit 1
  fi
fi

# Index output BAM or CRAM
if [[ -f "${output_prefix}.bam" ]]; then
  samtools index -@ "$nt" "${output_prefix}.bam" || exit 1

  # Check BAM EOF
  py_script="
import sys, os
def check_EOF(filename):
    EOF_hex = b'\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00'
    size = os.path.getsize(filename)
    with open(filename, 'rb') as f:
        f.seek(size - 28)
        if f.read(28) != EOF_hex:
            sys.stderr.write('EOF is missing\\n')
            sys.exit(1)
        else:
            sys.stderr.write('EOF is present\\n')
check_EOF('${output_prefix}.bam')
"
  python -c "$py_script" || exit 1
else
  mv "${input_file}" "${output_prefix}.cram" || exit 1
  samtools index -@ "$nt" "${output_prefix}.cram" || exit 1
fi

# Print output files
echo "[INFO] Finished. Output:"
[[ -f "${output_prefix}.bam" ]] && {
  echo "  - ${output_prefix}.bam"
  echo "  - ${output_prefix}.bam.bai"
}
[[ -f "${output_prefix}.cram" ]] && {
  echo "  - ${output_prefix}.cram"
  echo "  - ${output_prefix}.cram.crai"
}
