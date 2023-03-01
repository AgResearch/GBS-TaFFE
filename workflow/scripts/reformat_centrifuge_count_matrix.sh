#!/bin/bash

#!/bin/bash

# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# Utiltiy script to reformat a kraken formatted count matrix.
# Intended use is with centrifuge count matrixes generated using RE-RRS_GTDB snakemake worhflow.

# Input is a kraken2 formatted count matrix

# Returns 2 files:
#   1. Reformatted count matrix
#   2. #id sample key

usage() { echo "Usage $0 [-i <kraken.format.count.matrix.txt>] [-p <output/path>]" 1>&2; exit 1; }

get_arguments () {

  while getopts ":hi:p:" opt; do
    case $opt in
      h)
        usage ;;
      i)
        count_matrix="${OPTARG}" ;;
      p)
        out_path="${OPTARG}" ;;
      \?)
        echo "Unknown paramter: -${OPTARG}"
        usage
        ;;
    esac
  done
  
  if [ $OPTIND -eq 1 ]; then 
  usage 
  fi

  shift "$((OPTIND-1))"
}

check_arguments () {
  if [[ -z "${count_matrix}" ]]; then
    echo "Error: -i <kraken.format.count.matrix.txt> not given." 1>&2;
    exit 1
  fi
    if [[ -z "${out_path}" ]]; then
    echo "Error: -o <output/path> not given." 1>&2;
    exit 1
  fi

}

check_exsistence () {
  if [[ ! -f "${count_matrix}" ]]; then
    echo "Error: ${count_matrix} does not exist."
    exit 1
  fi

  if [[ ! -d "${out_path}" ]]; then
    echo "Error: ${out_path} does not exsist." 1>&2;
    exit 1
  fi

}

substitute_ids () {
  sample_id_key="${out_path}"/sample.id.key.txt
  if [[ -f "${sample_id_key}" ]]; then
    echo "Warning: Overwriting old ${sample_id_key}."
    echo > "${sample_id_key}"
  fi
  cat "${count_matrix}" | grep "#S" | sed "s/#S//g" | column -t | tr -s " " | awk '{gsub("[^:]*/", "", $2); print $0}' > "${sample_id_key}"
  
  updated_count_matrix="${out_path}"/clean.count.matrix.txt
  if [[ -f "${updated_count_matrix}" ]]; then
    echo "Warning: Overwriting old ${updated_count_matrix}."
    echo > "${updated_count_matrix}"
  fi
  printf 's/\t%s_/\t%s_/g\n' $(<"${sample_id_key}") > .patterns.sed
  cat "${count_matrix}" | grep "#perc" | sed "s/#//g" | sed -f .patterns.sed > "${updated_count_matrix}"
  cat "${count_matrix}" | grep -v "#" >> "${updated_count_matrix}"

}

main () {

  get_arguments "$@"
  check_arguments
  check_exsistence

  substitute_ids

}

main "$@"
