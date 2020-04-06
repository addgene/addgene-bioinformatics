from argparse import ArgumentParser
from configparser import ConfigParser
import csv
import os
import pickle

import utilities

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


def read_qc_sequences(sequencing_data_dir, qc_sequences_fNm):
    """
    Reads Addgene curated quality control full public sequences.

    Parameters
    ----------
    sequencing_data_dir : str
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
    data_pth = os.path.join(sequencing_data_dir, qc_sequences_fNm)
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


def align_cp_to_qc_sequences(assembler,
                             assembler_data_dir,
                             assembler_run_dir,
                             sequencing_data_dir,
                             qc_sequences_fNm,
                             aligner_config):
    """
    Aligns sequences produced by candidate sequence assembly processes
    to curated quality control full public sequences.

    Parameters
    ----------
    assembler : str
        name of assembler
    assembler_data_dir : str
        directory containing QC sequences file and candidate process
        run directories
    assembler_run_dir : str
        directory containing candidate process runs
    sequencing_data_dir : str
        directory containing QC sequences file
    qc_sequences_fNm : str
        file containing QC sequences
    aligner_config : dct
        pairwise aligner configuration

    Returns
    -------
    dct
        assembler, and apc if relevant, sequence, sequence length,
        doubled sequence, and corresponding alignment scores for each
        well in each plate in the run directory
    """
    # Create a pairwise aligner with the specified configuration
    aligner = utilities.create_aligner(aligner_config)

    # Read quality control (QC) sequences
    qc_sequences = read_qc_sequences(sequencing_data_dir, qc_sequences_fNm)

    # Initialize candidate process (CP) sequences dictionary
    cp_sequences = {}
    cp_sequences['aligner_config'] = aligner_config

    # Print header to a file, and stdout
    output_file = open(
        os.path.join(
            assembler_data_dir,
            assembler_run_dir + "-" + os.path.basename(aligner_config['file'])
            ).replace(".cfg", ".csv"), 'w')

    result = ""

    # TODO: Use more sensible formats
    result += "{:>11s},".format("assembler")
    result += "{:>8s},".format("plate")
    result += "{:>5s},".format("well")

    result += "{:>11s},".format("qc_seq_len")
    result += "{:>11s},".format("cp_seq_len")
    result += "{:>11s},".format("cp_seq_scr")
    result += "{:>11s},".format("cp_rcm_scr")

    result += "{:>12s},".format("apc_seq_len")
    result += "{:>12s},".format("apc_seq_scr")
    result += "{:>12s}".format("apc_rcm_scr")

    output_file.write(result + '\n')

    print(result)

    # Consider each plate directory contained in the run directory
    plate_specs = sorted(
        filter(
            '.DS_Store'.__ne__,
            os.listdir(
                os.path.join(assembler_data_dir, assembler_run_dir)
            )
        )
    )
    for plate_spec in plate_specs:

        # Find the plate id
        plate = plate_spec.split('_')[0]

        # Initialize plate dictionary
        if plate not in cp_sequences:
            cp_sequences[plate] = {}

        # Consider each well directory contained in the plate directory
        wells = sorted(
            filter(
                '.DS_Store'.__ne__,
                os.listdir(
                    os.path.join(
                        assembler_data_dir, assembler_run_dir, plate_spec
                    )
                )
            )
        )
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
            cp_sequences[plate][well]['qc'] = {}
            cp_sequences[plate][well][assembler] = {}
            if assembler == 'masurca':
                cp_sequence_fNm = "final.genome.scf.fasta"

            elif assembler == 'novoplasty':
                cp_sequence_fNm = "Circularized_assembly_1_Toil.fasta"

            elif assembler == 'shovill':
                cp_sequence_fNm = "contigs.fa"

            elif assembler == 'skesa':
                cp_sequence_fNm = "contigs.fa"

            elif assembler == 'spades':
                cp_sequence_fNm = "contigs.fasta"

            elif assembler == 'unicycler':
                cp_sequence_fNm = "assembly.fasta"

            elif assembler == 'seqwell':
                cp_sequence_fNm = "assembly.fasta"

            else:
                raise(Exception("Unexpected assembler"))
            try:
                cp_sequence = next(SeqIO.parse(
                    os.path.join(assembler_data_dir,
                                 assembler_run_dir,
                                 plate_spec,
                                 well,
                                 cp_sequence_fNm),
                    "fasta", IUPAC.unambiguous_dna)).seq
                cp_reverse_complement = cp_sequence.complement()[::-1]

                cp_sequence_score = aligner.score(
                    qc_doubled_sequence, cp_sequence)
                cp_reverse_complement_score = aligner.score(
                    qc_doubled_sequence, cp_reverse_complement)


            except Exception as e:
                cp_sequence = Seq("", IUPAC.unambiguous_dna)
                cp_reverse_complement = Seq("", IUPAC.unambiguous_dna)

                cp_sequence_score = 0
                cp_reverse_complement_score = 0

            cp_sequences[plate][well]['qc']['sequence'] = qc_sequence

            cp_sequences[plate][well][assembler]['sequence'] = cp_sequence
            cp_sequences[plate][well][assembler]['reverse_complement'] = cp_reverse_complement

            cp_sequences[plate][well][assembler]['sequence_score'] = cp_sequence_score
            cp_sequences[plate][well][assembler]['reverse_complement_score'] = cp_reverse_complement_score

            # Assign candidate process apc sequence values
            # Note: apc is not called after every assembler
            cp_sequences[plate][well]['apc'] = {}

            apc_sequence = Seq("", IUPAC.unambiguous_dna)
            apc_reverse_complement = Seq("", IUPAC.unambiguous_dna)

            apc_sequence_score = 0
            apc_reverse_complement_score = 0

            if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
                try:
                    apc_sequence = next(SeqIO.parse(
                        os.path.join(assembler_data_dir,
                                     assembler_run_dir,
                                     plate_spec,
                                     well,
                                     "apc.1.fa"),
                        "fasta", IUPAC.unambiguous_dna)).seq
                    apc_reverse_complement = apc_sequence.complement()[::-1]

                    apc_sequence_score = aligner.score(
                        qc_doubled_sequence, apc_sequence)
                    apc_reverse_complement_score = aligner.score(
                        qc_doubled_sequence, apc_reverse_complement)

                except Exception as e:
                    apc_sequence = Seq("", IUPAC.unambiguous_dna)
                    apc_reverse_complement = Seq("", IUPAC.unambiguous_dna)

                    apc_sequence_score = 0
                    apc_reverse_complement_score = 0

            cp_sequences[plate][well]['apc']['sequence'] = apc_sequence
            cp_sequences[plate][well]['apc']['reverse_complement'] = apc_reverse_complement

            cp_sequences[plate][well]['apc']['sequence_score'] = apc_sequence_score
            cp_sequences[plate][well]['apc']['reverse_complement_score'] = apc_reverse_complement_score

            # Print results to a file, and stdout
            result = ""

            result += "{:>11s},".format(assembler)
            result += "{:>8s},".format(plate)
            result += "{:>5s},".format(well)

            result += "{:11d},".format(qc_sequence_len)
            result += "{:11d},".format(len(cp_sequence))
            result += "{:11.1f},".format(cp_sequence_score)
            result += "{:11.1f},".format(cp_reverse_complement_score)

            result += "{:12d},".format(len(apc_sequence))
            result += "{:12.1f},".format(apc_sequence_score)
            result += "{:12.1f}".format(apc_reverse_complement_score)

            output_file.write(result + '\n')

            print(result)

    output_file.close()

    return cp_sequences


def accumulate_alignment_scores(assembler,
                                cp_sequences):
    """Accumulate alignement results.

    Finds the maximum of the sequence and reverse complement
    scores. Defines valid scores as maximum scores for which the
    assembler length is greater than zero. All results are converted
    to Numpy arrays.

    Parameters
    ----------
    assembler : str
        the assembler used to create candidate process sequences
    cp_sequences : dct
        assembler, and apc if relevant, sequence, sequence length,
        doubled sequence, and corresponding alignment scores for each
        well in each plate in the run directory

    Returns
    -------
    dct
         the input dictionary with sequence length, and random and
         maximum scores
    """
    # Initialize return values
    qc_sequence_len = []

    assembler_sequence_len = []
    assembler_maximum_score = []
    assembler_valid_index = []
    assembler_valid_score = []

    circularizer_sequence_len = []
    circularizer_maximum_score = []
    circularizer_valid_index = []
    circularizer_valid_score = []

    for plate, wells in cp_sequences.items():
        if plate == 'aligner_config' or not wells:
            continue

        for well, sequences in wells.items():
            if not sequences:
                continue

            qc_sequence_len.append(
                len(sequences['qc']['sequence']))

            assembler_sequence_len.append(
                len(sequences[assembler]['sequence']))
            assembler_maximum_score.append(max(
                [sequences[assembler]['sequence_score'],
                 sequences[assembler]['reverse_complement_score']]))

            if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
                circularizer_sequence_len.append(
                    len(sequences['apc']['sequence']))

                circularizer_maximum_score.append(max(
                    [sequences['apc']['sequence_score'],
                     sequences['apc']['reverse_complement_score']]))

    qc_sequence_len = np.array(qc_sequence_len)

    assembler_sequence_len = np.array(assembler_sequence_len)
    assembler_maximum_score = np.array(assembler_maximum_score)
    assembler_valid_score = np.zeros(assembler_maximum_score.shape)

    circularizer_sequence_len = np.array(circularizer_sequence_len)
    circularizer_maximum_score = np.array(circularizer_maximum_score)
    circularizer_valid_score = np.zeros(assembler_maximum_score.shape)

    # Identify results for which the sequence length is greater than
    # zero
    assembler_valid_index = (assembler_sequence_len > 0)
    assembler_valid_score[assembler_valid_index] = (
        assembler_maximum_score[assembler_valid_index])
    if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
        circularizer_valid_index = (circularizer_sequence_len > 0)
        circularizer_valid_score[circularizer_valid_index] = (
            circularizer_maximum_score[circularizer_valid_index])

    # Assign return values
    cp_sequences['qc_sequence_len'] = qc_sequence_len

    cp_sequences['assembler_sequence_len'] = assembler_sequence_len
    cp_sequences['assembler_maximum_score'] = assembler_maximum_score
    cp_sequences['assembler_valid_index'] = assembler_valid_index
    cp_sequences['assembler_valid_score'] = assembler_valid_score

    cp_sequences['circularizer_sequence_len'] = circularizer_sequence_len
    cp_sequences['circularizer_maximum_score'] = circularizer_maximum_score
    cp_sequences['circularizer_valid_index'] = circularizer_valid_index
    cp_sequences['circularizer_valid_score'] = circularizer_valid_score

    return cp_sequences


def plot_alignment_scores(assembler, cp_sequences):
    """Plot histograms of relative alignment scores resulting form
    assembled and circularized sequences.

    Parameters
    ----------
    assembler : str
        the assembler used to create candidate process sequences
    cp_sequences : dct
        assembler, and apc if relevant, sequence, sequence length,
        doubled sequence, and corresponding alignment scores for each
        well in each plate in the run directory
    """
    fig, ax = plt.subplots()

    bin_width = 5.0
    bin_edges = np.arange(0.0, 100.0 + 2 * bin_width, bin_width)

    # Compute and plot assembler score
    assembler_valid_index = cp_sequences['assembler_valid_index']
    assembler_valid_score = (
        100.0 * (1 +
                 (cp_sequences['assembler_valid_score'][assembler_valid_index]
                  - cp_sequences['qc_sequence_len'][assembler_valid_index])
                 / cp_sequences['qc_sequence_len'][assembler_valid_index])
   )

    assembler_histogram, tmp = np.histogram(assembler_valid_score, bin_edges)
    ax.bar(bin_edges[:-1] - bin_width / 4,
           100.0 * assembler_histogram / len(assembler_valid_index),
           width=2.0, align='center', color='b')

    # Compute and plot circularizer score
    if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
        circularizer_valid_index = cp_sequences['circularizer_valid_index']
        circularizer_valid_score = (
            100.0 * (1 +
                     (cp_sequences['circularizer_valid_score'][circularizer_valid_index]
                      - cp_sequences['qc_sequence_len'][circularizer_valid_index])
                     / cp_sequences['qc_sequence_len'][circularizer_valid_index])
        )

        circularizer_histogram, tmp = np.histogram(circularizer_valid_score, bin_edges)
        ax.bar(bin_edges[:-1] + bin_width / 4,
               100.0 * circularizer_histogram / len(circularizer_valid_index),
               width=2.0, align='center', color='r')

    ax.set_title(assembler)
    ax.set_xlabel("Alignment Score (% of QC Sequence Length)")
    ax.set_ylabel("Fraction of Assemblies (%)")

    x1, x2 = ax.get_xlim()
    xW = x2 - x1

    y1, y2 = ax.get_ylim()
    yW = y2 - y1

    ax.text(x1 + 0.10 * xW, y2 - 0.05 * yW,
            "Assembled: {}".format(len(assembler_valid_index)))

    plt.show()


if __name__ == "__main__":

    # Add and parse arguments
    parser = ArgumentParser()
    base_dir = os.path.join(
        os.sep, *os.path.abspath(__file__).split(os.sep)[0:-4])
    home_dir = os.path.join(
        os.sep, *os.path.abspath(__file__).split(os.sep)[0:-3])
    parser.add_argument(
        '-d', '--assembler-data-directory',
        default=os.path.join(
            base_dir,
            "addgene-assembler-data",
            "results-2020-01-16"),
        help="directory containing assembler run directory")
    parser.add_argument(
        '-s', '--sequencing-data-directory',
        default=os.path.join(
            base_dir,
            "addgene-sequencing-data"),
        help="directory containing sequencing data directories")
    parser.add_argument(
        '-q', '--qc-sequences-file',
        default="2018_QC_Sequences.csv",
        help="file containing QC sequences")
    parser.add_argument(
        '-p', '--plot',
        action='store_true',
        help="plot alignment results")
    parser.add_argument(
        '-r', '--reprocess',
        action='store_true',
        help="reprocess each run directory")
    options = parser.parse_args()

    # Read and parse configuration file
    config_dir = os.path.join(home_dir, "resources")
    config_file = "pairwise-rl-01.cfg"
    config = ConfigParser()
    config.read(os.path.join(config_dir, config_file))
    config['aligner']['file'] = config_file

    # Identify run directories
    assembler_run_dirs = []
    for assembler_run_dir in os.listdir(options.assembler_data_directory):
        if os.path.isdir(
            os.path.join(
                options.assembler_data_directory,
                assembler_run_dir)):
            assembler_run_dirs.append(assembler_run_dir)

    # Initialize candidate process alignment dictionary
    cp_alignment = {}
    cp_alignment['aligner'] = dict(config['aligner'])
    for assembler_run_dir in assembler_run_dirs:

        # Identify assembler corresponding to the run directory
        assembler = assembler_run_dir.split('-')[0]

        # Initialize candidate process assembler dictionary
        cp_alignment[assembler] = {}

        # Create and dump, or load candidate process sequences pickle
        cp_alignment_path = os.path.join(
            options.assembler_data_directory,
            assembler_run_dir + "-" + os.path.basename(config_file)
            ).replace(".cfg", ".pickle")
        if not os.path.isfile(cp_alignment_path) or options.reprocess:

            # Align sequences produced by candidate sequence assembly
            # processes to curated quality control full public
            # sequences
            cp_sequences = align_cp_to_qc_sequences(
                assembler,
                options.assembler_data_directory,
                assembler_run_dir,
                options.sequencing_data_directory,
                options.qc_sequences_file,
                dict(config['aligner']))

            with open(cp_alignment_path, 'wb') as f:
                pickle.dump(cp_sequences, f, pickle.HIGHEST_PROTOCOL)

        else:

            with open(cp_alignment_path, 'rb') as f:
                cp_sequences = pickle.load(f)

        # Accumulate alignement results
        cp_sequences = accumulate_alignment_scores(assembler, cp_sequences)
        cp_alignment[assembler]['sequences'] = cp_sequences

        # Note in computing valid indexes below that since the
        # assembler results are from the same plates, the nonzero QC
        # sequences are the same for each assembler
        qc_valid_index = (cp_sequences['qc_sequence_len'] > 0)

        # TODO: Review all following lines

        # Plot alignment results
        if options.plot:
            plot_alignment_scores(assembler, cp_sequences)

        # Print alignment results
        print("")
        print("Assembler: {}".format(assembler))
        print("    Assembled (nonzero length): {:.1f}%".format(
            100 * sum(cp_sequences['assembler_sequence_len'] > 0)
            / len(cp_sequences['assembler_sequence_len'])))
        if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
            print("    Circularized (nonzero length): {:.1f}%".format(
                100 * sum(cp_sequences['circularizer_sequence_len'] > 0)
                / sum(cp_sequences['assembler_sequence_len'] > 0)))

        fraction_aligned = 0.10
        print("    Assembly aligned (within {:.1f}%): {:.1f}%".format(
            100 * fraction_aligned,
            100 * sum(
                (cp_sequences['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                 cp_sequences['assembler_valid_score'][qc_valid_index]) &
                (cp_sequences['assembler_valid_score'][qc_valid_index] <=
                 cp_sequences['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned)))
            / len(cp_sequences['assembler_sequence_len'][qc_valid_index])))
        if assembler in utilities.ASSEMBLERS_REQUIRING_APC:
            print("    Circularized assembly aligned (within {:.1f}%): {:.1f}%".format(
                100 * fraction_aligned,
                100 * sum(
                    (cp_sequences['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                     cp_sequences['circularizer_valid_score'][qc_valid_index]) &
                    (cp_sequences['circularizer_valid_score'][qc_valid_index] <=
                     cp_sequences['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned)))
                / sum(cp_sequences['assembler_sequence_len'][qc_valid_index] > 0)))

    if assembler != "seqwell":
        print("")
        print("Assembler: any")
        print("    Assembled (nonzero length): {:.1f}%".format(
            100 * sum(
                # (cp_alignment['masurca']['sequences']['assembler_sequence_len'] > 0) |
                (cp_alignment['novoplasty']['sequences']['assembler_sequence_len'] > 0) |
                (cp_alignment['shovill']['sequences']['assembler_sequence_len'] > 0) |
                # (cp_alignment['skesa']['sequences']['assembler_sequence_len'] > 0) |
                (cp_alignment['spades']['sequences']['assembler_sequence_len'] > 0) # |
                # (cp_alignment['unicycler']['sequences']['assembler_sequence_len'] > 0)
            ) / len(cp_sequences['assembler_sequence_len'])))
        print("    Circularized (nonzero length): {:.1f}%".format(
            100 * sum(
                # (cp_alignment['masurca']['sequences']['circularizer_sequence_len'] > 0) |
                (cp_alignment['novoplasty']['sequences']['assembler_sequence_len'] > 0) |
                (cp_alignment['shovill']['sequences']['circularizer_sequence_len'] > 0) |
                # (cp_alignment['skesa']['sequences']['circularizer_sequence_len'] > 0) |
                (cp_alignment['spades']['sequences']['circularizer_sequence_len'] > 0) # |
                # (cp_alignment['unicycler']['sequences']['circularizer_sequence_len'] > 0)
            ) / len(cp_sequences['assembler_sequence_len'])))

        print("    Linear assembly aligned (within {:.1f}%): {:.1f}".format(
            100 * fraction_aligned,
            100 * sum(
                # ((cp_alignment['masurca']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['masurca']['sequences']['assembler_valid_score'][qc_valid_index]) &
                #  (cp_alignment['masurca']['sequences']['assembler_valid_score'][qc_valid_index] <=
                #   cp_alignment['masurca']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['novoplasty']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['novoplasty']['sequences']['assembler_valid_score'][qc_valid_index]) &
                 (cp_alignment['novoplasty']['sequences']['assembler_valid_score'][qc_valid_index] <=
                  cp_alignment['novoplasty']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['shovill']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['shovill']['sequences']['assembler_valid_score'][qc_valid_index]) &
                 (cp_alignment['shovill']['sequences']['assembler_valid_score'][qc_valid_index] <=
                  cp_alignment['shovill']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                # ((cp_alignment['skesa']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['skesa']['sequences']['assembler_valid_score'][qc_valid_index]) &
                #  (cp_alignment['skesa']['sequences']['assembler_valid_score'][qc_valid_index] <=
                #   cp_alignment['skesa']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['spades']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['spades']['sequences']['assembler_valid_score'][qc_valid_index]) &
                 (cp_alignment['spades']['sequences']['assembler_valid_score'][qc_valid_index] <=
                  cp_alignment['spades']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) # |

                # ((cp_alignment['unicycler']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['unicycler']['sequences']['assembler_valid_score'][qc_valid_index]) &
                #  (cp_alignment['unicycler']['sequences']['assembler_valid_score'][qc_valid_index] <=
                #   cp_alignment['unicycler']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned)))

            ) / len(cp_alignment['novoplasty']['sequences']['assembler_sequence_len'])))

        print("    Circular assembly aligned (within {:.1f}%): {:.1f}".format(
            100 * fraction_aligned,
            100 * sum(
                # ((cp_alignment['masurca']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['masurca']['sequences']['circularizer_valid_score'][qc_valid_index]) &
                #  (cp_alignment['masurca']['sequences']['circularizer_valid_score'][qc_valid_index] <=
                #   cp_alignment['masurca']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['novoplasty']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['novoplasty']['sequences']['assembler_valid_score'][qc_valid_index]) &
                 (cp_alignment['novoplasty']['sequences']['assembler_valid_score'][qc_valid_index] <=
                  cp_alignment['novoplasty']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['shovill']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['shovill']['sequences']['circularizer_valid_score'][qc_valid_index]) &
                 (cp_alignment['shovill']['sequences']['circularizer_valid_score'][qc_valid_index] <=
                  cp_alignment['shovill']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                # ((cp_alignment['skesa']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['skesa']['sequences']['circularizer_valid_score'][qc_valid_index]) &
                #  (cp_alignment['skesa']['sequences']['circularizer_valid_score'][qc_valid_index] <=
                #   cp_alignment['skesa']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) |

                ((cp_alignment['spades']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                  cp_alignment['spades']['sequences']['circularizer_valid_score'][qc_valid_index]) &
                 (cp_alignment['spades']['sequences']['circularizer_valid_score'][qc_valid_index] <=
                  cp_alignment['spades']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned))) # |

                # ((cp_alignment['unicycler']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 - fraction_aligned) <=
                #   cp_alignment['unicycler']['sequences']['circularizer_valid_score'][qc_valid_index]) &
                #  (cp_alignment['unicycler']['sequences']['circularizer_valid_score'][qc_valid_index] <=
                #   cp_alignment['unicycler']['sequences']['qc_sequence_len'][qc_valid_index] * (1.0 + fraction_aligned)))

            ) / len(cp_alignment['novoplasty']['sequences']['assembler_sequence_len'][qc_valid_index])))
