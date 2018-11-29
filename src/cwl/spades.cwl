#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: A single SPAdes run

doc: |
  Accepts paired-end Illumina reads for assembly with a coverage
  cutoff specification

baseCommand: spades.py

inputs:

  read_one_file:
    type: File
    label: FASTQ Illumina short left paired reads
    format: edam:format_1931
    inputBinding:
      prefix: "-1"
      separate: true
      position: 1

  read_two_file:
    type: File
    label: FASTQ Illumina short right paired reads
    format: edam:format_1931
    inputBinding:
      prefix: "-2"
      separate: true
      position: 2

  output_directory:
    type: string
    inputBinding:
      position: 3
      prefix: -o
      separate: true

  coverage_cutoff:
    type: int
    inputBinding:
      position: 4
      prefix: --cov-cutoff
      separate: true

outputs:

  log:
    type: File
    label: Run log (same as stdout)
    format: edam:format_3671
    outputBinding:
      glob: $(inputs.output_directory)/spades.log

  warnings:
    type: File
    label: Run warnings
    format: edam:format_3671
    outputBinding:
      glob: $(inputs.output_directory)/warnings.log

  scaffolds:
    type: File
    label: Assembled scaffolds
    format: edam:format_1929
    outputBinding:
      glob: $(inputs.output_directory)/scaffolds.fasta

$namespaces:
  edam: http://edamontology.org/

$schemas:
  - http://edamontology.org/EDAM_1.18.owl
