#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: A single apc run

doc: |
  Tests each sequence in a supplied fasta file for self-overlap. If
  overlap is found, the 5' copy of the overlapping sequence is
  trimmed, the ends joined, and the join moved to the middle of the
  output sequence. The join location is appended to the sequence
  header, to allow confirmation.

baseCommand: apc.pl

stdout: apc.out

inputs:

  basename:
    type: string
    label: Basename for output files
    inputBinding:
      position: 1
      prefix: -b
      separate: true

  scaffolds:
    type: File
    label: Assembled scaffolds
    format: edam:format_1929
    inputBinding:
      position: 2

outputs:

  output:
    type: stdout

  alignments:
    type:
      type: array
      items: File
    label: LAST output
    format: edam:format_3671
    outputBinding:
      glob: $(inputs.basename)_aln_*.txt

  sequence:
    type: File
    label: Circularized sequence
    format: edam:format_1929
    outputBinding:
      glob: $(inputs.basename).1.fa

$namespaces:
  edam: http://edamontology.org/

$schemas:
  - http://edamontology.org/EDAM_1.18.owl
