from argparse import ArgumentParser
import csv
import os
import pickle

from Bio import Align
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from matplotlib import pyplot as plt
import numpy as np


def to_int(a_str):
    """
    Convert a string to an integer, return 0 for empty strings.

    Parameters
    ----------
    a_str : str
        a str

    Returns
    -------
    int
        the int
    """
    if a_str == '':
        return 0
    else:
        return int(a_str)


def read_qc_sequences(assembler_data_dir, qc_sequences_fNm):
    """
    Reads Addgene curated quality control full public sequences.

    Parameters
    ----------
    assembler_data_dir : str
        directory containing QC sequences file
    qc_sequences_fNm : str
        file containing QC sequences

    Returns
    -------
    dict
        plasmid and sequences ids, sequence, sequence length, and
        doubled sequence
    """
    # Read QC sequences
    qc_sequences = {}
    data_pth = os.path.join(assembler_data_dir, qc_sequences_fNm)
    with open(data_pth) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:

            # Sequencing Plate barcode
            plate = row[1]
            if plate not in qc_sequences:
                qc_sequences[plate] = {}

            # Well (location in 96-well plate)
            well = row[2]
            if well not in qc_sequences[plate]:
                qc_sequences[plate][well] = {}

            # Plasmid ID
            qc_sequences[plate][well]['plasmid_id'] = to_int(row[0])

            # Sequence ID (LIMS identifier)
            qc_sequences[plate][well]['sequence_id'] = to_int(row[3])

            # Sequence length (of QC-team approved sequence)
            qc_sequences[plate][well]['sequence_len'] = to_int(row[4])

            # Sequence (actual nucleotide sequence)
            qc_sequences[plate][well]['sequence'] = Seq(row[5])

            # Doubled sequence (actual nucleotide sequence, but doubled)
            qc_sequences[plate][well]['doubled_sequence'] = Seq(row[5]
                                                                + row[5])

    return qc_sequences


def align_cp_to_qc_sequences(assembler_data_dir,
                             qc_sequences_fNm,
                             assembler_run_dir):
    """
    Aligns sequences produced by candidate sequence assembly processes
    to curated quality control full public sequences.

    Parameters
    ----------
    assembler_data_dir : str
        directory containing QC sequences file and candidate process
        run directories
    qc_sequences_fNm : str
        file containing QC sequences
    assembler_run_dir : str
        directory containing candidate process runs

    Returns
    -------
    dict
        assembler, and apc if relevant, sequence, sequence length,
        doubled sequence, and corresponding alignment scores for each
        well in each plate in the run directory
    """

    # Identify assembler corresponding to the run directory
    assembler = assembler_run_dir.split('-')[0]

    # Create a pairwise aligner with these parameters:
    #     algorithm: 'Needleman-Wunsch'
    #     match_score: 1.000000
    #     mismatch_score: 0.000000
    #     target_open_gap_score: 0.000000
    #     target_extend_gap_score: 0.000000
    #     target_left_open_gap_score: 0.000000
    #     target_left_extend_gap_score: 0.000000
    #     target_right_open_gap_score: 0.000000
    #     target_right_extend_gap_score: 0.000000
    #     query_open_gap_score: 0.000000
    #     query_extend_gap_score: 0.000000
    #     query_left_open_gap_score: 0.000000
    #     query_left_extend_gap_score: 0.000000
    #     query_right_open_gap_score: 0.000000
    #     query_right_extend_gap_score: 0.000000
    #     mode: global
    aligner = Align.PairwiseAligner()

    # Read quality control (QC) sequences
    qc_sequences = read_qc_sequences(assembler_data_dir, qc_sequences_fNm)

    # Initialize candidate process (CP) sequences dictionary
    cp_sequences = {}

    # Consider each plate directory contained in the run directory
    output_file = open(
        os.path.join(assembler_data_dir, assembler_run_dir + ".csv"), 'w')

    plate_specs = sorted(
        os.listdir(
            os.path.join(assembler_data_dir, assembler_run_dir)))
    for plate_spec in plate_specs:

        # Find the plate id
        plate = plate_spec.split('_')[0]

        # Initialize plate dictionary
        if plate not in cp_sequences:
            cp_sequences[plate] = {}

        # Consider each well directory contained in the plate directory
        wells = sorted(
            os.listdir(
                os.path.join(assembler_data_dir, assembler_run_dir, plate_spec)))
        for well in wells:

            # Initialize well dictionary
            if well not in cp_sequences[plate]:
                cp_sequences[plate][well] = {}

            # Assign quality control sequence values
            try:
                qc_sequence = qc_sequences[plate][well]['sequence']
                qc_sequence_len = qc_sequences[plate][well]['sequence_len']
                qc_doubled_sequence = qc_sequences[plate][well]['doubled_sequence']
            except Exception as e:
                continue

            # Assign candidate process assembler sequence values
            cp_sequences[plate][well][assembler] = {}
            cp_sequence_fNm = ""
            if assembler == "novoplasty":
                cp_sequence_fNm = "Circularized_assembly_1_Toil.fasta"

            elif assembler == "shovill":
                cp_sequence_fNm = "contigs.fa"

            elif assembler == "spades":
                cp_sequence_fNm = "contigs.fasta"
            try:
                cp_sequence = next(SeqIO.parse(
                    os.path.join(assembler_data_dir,
                                 assembler_run_dir,
                                 plate_spec,
                                 well,
                                 cp_sequence_fNm),
                    "fasta", IUPAC.unambiguous_dna)).seq
                cp_complement = cp_sequence.complement()
                cp_reverse = cp_sequence[len(cp_sequence):0:-1]
                cp_reverse_complement = cp_complement[len(cp_complement):0:-1]

                cp_sequence_score = aligner.score(
                    qc_doubled_sequence, cp_sequence)
                cp_complement_score = aligner.score(
                    qc_doubled_sequence, cp_complement)
                cp_reverse_score = aligner.score(
                    qc_doubled_sequence, cp_reverse)
                cp_reverse_complement_score = aligner.score(
                    qc_doubled_sequence, cp_reverse_complement)

            except Exception as e:
                cp_sequence = Seq("", IUPAC.unambiguous_dna)
                cp_complement = Seq("", IUPAC.unambiguous_dna)
                cp_reverse = Seq("", IUPAC.unambiguous_dna)
                cp_reverse_complement = Seq("", IUPAC.unambiguous_dna)

                cp_sequence_score = 0
                cp_complement_score = 0
                cp_reverse_score = 0
                cp_reverse_complement_score = 0

            cp_sequences[plate][well][assembler]['sequence'] = cp_sequence
            cp_sequences[plate][well][assembler]['complement'] = cp_complement
            cp_sequences[plate][well][assembler]['reverse'] = cp_reverse
            cp_sequences[plate][well][assembler]['reverse_complement'] = cp_reverse_complement

            cp_sequences[plate][well][assembler]['sequence_score'] = cp_sequence_score
            cp_sequences[plate][well][assembler]['complement_score'] = cp_complement_score
            cp_sequences[plate][well][assembler]['reverse_score'] = cp_reverse_score
            cp_sequences[plate][well][assembler]['reverse_complement_score'] = cp_reverse_complement_score

            # Assign candidate process apc sequence values
            # Note: apc is not called after every assembler
            cp_sequences[plate][well]['apc'] = {}

            apc_sequence = Seq("", IUPAC.unambiguous_dna)
            apc_complement = Seq("", IUPAC.unambiguous_dna)
            apc_reverse = Seq("", IUPAC.unambiguous_dna)
            apc_reverse_complement = Seq("", IUPAC.unambiguous_dna)

            apc_sequence_score = 0
            apc_complement_score = 0
            apc_reverse_score = 0
            apc_reverse_complement_score = 0

            if assembler in ["shovill", "spades"]:
                try:
                    apc_sequence = next(SeqIO.parse(
                        os.path.join(assembler_data_dir,
                                     assembler_run_dir,
                                     plate_spec,
                                     well,
                                     "apc.1.fa"),
                        "fasta", IUPAC.unambiguous_dna)).seq
                    apc_complement = apc_sequence.complement()
                    apc_reverse = apc_sequence[len(apc_sequence):0:-1]
                    apc_reverse_complement = apc_complement[len(apc_complement):0:-1]

                    apc_sequence_score = aligner.score(
                        qc_doubled_sequence, apc_sequence)
                    apc_complement_score = aligner.score(
                        qc_doubled_sequence, apc_complement)
                    apc_reverse_score = aligner.score(
                        qc_doubled_sequence, apc_reverse)
                    apc_reverse_complement_score = aligner.score(
                        qc_doubled_sequence, apc_reverse_complement)

                except Exception as e:
                    apc_sequence = Seq("", IUPAC.unambiguous_dna)
                    apc_complement = Seq("", IUPAC.unambiguous_dna)
                    apc_reverse = Seq("", IUPAC.unambiguous_dna)
                    apc_reverse_complement = Seq("", IUPAC.unambiguous_dna)

                    apc_sequence_score = 0
                    apc_complement_score = 0
                    apc_reverse_score = 0
                    apc_reverse_complement_score = 0

            cp_sequences[plate][well]['apc']['sequence'] = apc_sequence
            cp_sequences[plate][well]['apc']['complement'] = apc_complement
            cp_sequences[plate][well]['apc']['reverse'] = apc_reverse
            cp_sequences[plate][well]['apc']['reverse_complement'] = apc_reverse_complement

            cp_sequences[plate][well]['apc']['sequence_score'] = apc_sequence_score
            cp_sequences[plate][well]['apc']['complement_score'] = apc_complement_score
            cp_sequences[plate][well]['apc']['reverse_score'] = apc_reverse_score
            cp_sequences[plate][well]['apc']['reverse_complement_score'] = apc_reverse_complement_score

            # Print results to a file, and stdout
            result = ""

            result += "{:10s}".format(assembler)
            result += ", {:7s}".format(plate)
            result += ", {:3s}".format(well)

            result += ", {:10s}".format(
                str(qc_sequence[0:min(10, len(qc_sequence))]))
            result += ", {:10s}".format(
                str(cp_sequence[0:min(10, len(cp_sequence))]))
            result += ", {:10s}".format(
                str(apc_sequence[0:min(10, len(apc_sequence))]))

            if len(qc_sequence) > 0:
                result += ", {:5.1f}".format(
                    100 * cp_sequence_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * cp_complement_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * cp_reverse_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * cp_reverse_complement_score / len(qc_sequence))

                result += ", {:5.1f}".format(
                    100 * apc_sequence_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * apc_complement_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * apc_reverse_score / len(qc_sequence))
                result += ", {:5.1f}".format(
                    100 * apc_reverse_complement_score / len(qc_sequence))

            else:
                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)

                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)
                result += ", {:5.1f}".format(0)

            output_file.write(result + '\n')
            print(result)

    output_file.close()

    return cp_sequences


if __name__ == "__main__":

    # And and parse arguments
    parser = ArgumentParser()
    parser.add_argument('-r', '--reprocess',
                        action='store_true',
                        help="reprocess each run directory")
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
    parser.add_argument('-d', '--assembler-data-directory',
                        default=os.path.join(os.sep + os.path.join(*cmps),
                                             "addgene-assembler-data"),
                        help="directory containing assembler run directory")
    parser.add_argument('-q', '--qc-sequences-file',
                        default=os.path.join(os.sep + os.path.join(*cmps),
                                             "addgene-assembler-data",
                                             "2018_QC_Sequences.csv"),
                        help="file containing QC sequences")
    options = parser.parse_args()

    # Identify run directories
    assembler_run_dirs = []
    for assembler_run_dir in os.listdir(options.assembler_data_directory):
        if os.path.isdir(
            os.path.join(options.assembler_data_directory,
                         assembler_run_dir)):
            assembler_run_dirs.append(assembler_run_dir)

    # Initialize candidate process alignment dictionary
    cp_alignment = {}

    for assembler_run_dir in assembler_run_dirs:
        assembler = assembler_run_dir.split('-')[0]

        # Initialize candidate process assembler dictionary
        cp_alignment[assembler] = {}

        # Create and dump, or load candidate process sequences pickle
        cp_alignment_pNm = assembler_run_dir + ".pickle"
        cp_alignment_pth = os.path.join(
            options.assembler_data_directory, cp_alignment_pNm)
        if not os.path.isfile(cp_alignment_pth) or options.reprocess:

            # Align sequences produced by candidate sequence assembly
            # processes to curated quality control full public
            # sequences
            cp_sequences = align_cp_to_qc_sequences(
                options.assembler_data_directory,
                options.qc_sequences_file,
                assembler_run_dir)

            with open(cp_alignment_pth, 'wb') as f:
                pickle.dump(cp_sequences, f, pickle.HIGHEST_PROTOCOL)

        else:

            with open(cp_alignment_pth, 'rb') as f:
                cp_sequences = pickle.load(f)

        cp_alignment[assembler] = cp_sequences

        """
        assembler_sequences = []
        circularizer_sequences = []

        for plate, wells in cp_sequences.items():
            if not wells:
                continue
            
            for well, sequences in wells.items():
                if not sequences:
                    continue
                assembler_sequences.append(sequences[assembler])
                if assembler in ["shovill", "spades"]:
                    circularizer_sequences.append(sequences['apc'])

    fig, axs = plt.subplots()
    """