#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow

label: A single SPAdes, then apc, run

doc: |
  SPAdes accepts paired-end Illumina reads for assembly with a
  coverage cutoff specification

  apc tests each sequence in a supplied fasta file for
  self-overlap. If overlap is found, the 5' copy of the overlapping
  sequence is trimmed, the ends joined, and the join moved to the
  middle of the output sequence. The join location is appended to the
  sequence header, to allow confirmation.

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

  basename:
    type: string
    label: Basename for output files
    inputBinding:
      position: 1
      prefix: -b
      separate: true

outputs:

  log:
    type: File
    label: Run log (same as stdout)
    format: edam:format_3671
    outputSource: spades/log

  warnings:
    type: File
    label: Run warnings
    format: edam:format_3671
    outputSource: spades/warnings

  scaffolds:
    type: File
    label: Assembled scaffolds
    format: edam:format_1929
    outputSource: spades/scaffolds

  output:
    type: File
    label: Run log (same as stdout)
    format: edam:format_3671
    outputSource: apc/output

  alignments:
    type: File[]
    label: LAST output
    format: edam:format_3671
    outputSource: apc/alignments

  sequence:
    type: File
    label: Circularized sequence
    format: edam:format_1929
    outputSource: apc/sequence
    
steps:

  spades:

    run: spades.cwl

    in:
      read_one_file: read_one_file
      read_two_file: read_two_file
      output_directory: output_directory
      coverage_cutoff: coverage_cutoff

    out: [log, warnings, scaffolds]

  apc:

    run: apc.cwl

    in:
      basename: basename
      scaffolds: spades/scaffolds

    out: [output, alignments, sequence]

$namespaces:
  edam: http://edamontology.org/

$schemas:
  - http://edamontology.org/EDAM_1.18.owl
