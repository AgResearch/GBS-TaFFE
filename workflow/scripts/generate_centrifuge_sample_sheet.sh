#!/bin/bash

# 2022 Benjamin J Perry
# MIT License
# Copyright (c) 2022 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

# Utiltiy script to generate sample sheet for centrifuge-multi.

# Input is a path to directory of fastq data, 
# a specific pattern to match with, 
# and a path to an output file to write to.

# Output is a 5 column tab-delimited file without a header line, in the format:

# Column 1 - [1/2] indicating single or paired reads
# Column 2 - read1 path
# Column 3 - read2 path [disregarded in single-end mode]
# Column 4 - centrifuge output file
# Column 5 - centrifuge report file



usage() { echo "Usage $0 [-d <results/directory>] [-p <file.name.suffix.fastq>] [-o <output/samplesheet.tsv>]" 1>&2; exit 1; }

get_arguments () {


  while getopts ":hd:p:o:" opt; do
    case $opt in
      h)
        usage ;;
      d)
        results_directory="${OPTARG}" ;;
      p)
        file_pattern="${OPTARG}" ;;
      o)
        sample_sheet="${OPTARG}" ;;
      \?)
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
  if [[ -z "${results_directory}" ]]; then
    echo "Error: -d <results/directory> not given." 1>&2;
    exit 1
  fi

  if [[ -z "${file_pattern}" ]]; then
    echo "Error: -p <file.name.suffix.fastq> not given." 1>&2;
    exit 1
  fi

  if [[ -z "${sample_sheet}" ]]; then
    echo "Error: -o <output/samplesheet.tsv> not given." 1>&2;
    exit 1
  fi
}

check_exsistence () {
  if [[ ! -d "${results_directory}" ]]; then
    echo "Error: ${results_directory} does not exist."
    exit 1
  fi

  matches="$(ls ${results_directory}/*${file_pattern})"
  if [[ -z "${matches}" ]]; then
    echo "Error: No files in ${results_directory} which match the pattern ${file_pattern}"
    exit 1
  fi
  
  if [[ -f "${sample_sheet}" ]]; then
    echo "Warning: ${sample_sheet} already exsists. Overwriting."
    echo > "${sample_sheet}"
  else
    mkdir -p $(dirname "${sample_sheet}") && touch "${sample_sheet}"
  fi

}

generate_centrifuge_sample_sheet () {
  for file in $matches; do
    fid="$(echo ${file} | xargs basename | cut -d "_" -f 1)"
    printf "1\t${file}\tNA\tresults/03_centrifuge/${fid}.GTDB.centrifuge\tresults/03_centrifuge/${fid}.GTDB.centrifuge.report\n" >> "${sample_sheet}"
    done

  echo "Completed centrifuge sample sheet."
  echo "${sample_sheet}: "
  cat "${sample_sheet}" | column -t | head
  echo "[...]"

}

main () {

    get_arguments "$@"

    check_arguments

    check_exsistence

    generate_centrifuge_sample_sheet

    exit 0
}


main "$@"
