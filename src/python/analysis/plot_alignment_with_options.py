from argparse import ArgumentParser
from configparser import ConfigParser
import os
import pickle
from random import choice, choices

import utilities

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import numpy as np


def create_o_seq(a_seq, n_off):
    """Offsets a sequence by the specified number of nucleotides.

    Parameters
    ----------
    a_seq : Bio.Seq.Seq
        the sequence to offset
    n_off : int
        the number of nucleotides offset

    Returns
    -------
    Bio.Seq.Seq
        the offset sequence
    """
    o_seq = a_seq[n_off:] + a_seq[:n_off]
    return o_seq


def create_e_seq(o_seq, n_err):
    """Replaces the specified number of nucleotides by a randomly
    chosen other nucleotide in the specified sequence.

    Parameters
    ----------
    o_seq : Bio.Seq.Seq
        the sequence in which to replace nucleotides
    n_err : int
        the number of nucleotides to replace

    Returns
    -------
    Bio.Seq.Seq
        the input sequence with errors
    """
    seq_len = len(o_seq)
    e_seq_lst = list(o_seq)
    for nt_idx in choices(range(seq_len), k=n_err):
        nt_lst = list("CGTA")
        del nt_lst[nt_lst.index(e_seq_lst[nt_idx])]
        e_seq_lst[nt_idx] = choice(nt_lst)
    e_seq_str = "".join(e_seq_lst)
    e_seq = Seq(e_seq_str, IUPAC.unambiguous_dna)
    return e_seq


def test_aligners():
    """Compare the scores produced by aligners corresponding to the
    default, and two recommended sets of parameters.

    """
    # Read and parse configuration files
    config_deflt = ConfigParser()
    config_deflt.read("../../resources/pairwise-default.cfg")
    config_jn_01 = ConfigParser()
    config_jn_01.read("../../resources/pairwise-jn-01.cfg")
    config_rl_01 = ConfigParser()
    config_rl_01.read("../../resources/pairwise-rl-01.cfg")

    # Create aligners using the specified parameters
    aligner_deflt = utilities.create_aligner(config_deflt["aligner"])
    aligner_jn_01 = utilities.create_aligner(config_jn_01["aligner"])
    aligner_rl_01 = utilities.create_aligner(config_rl_01["aligner"])

    # Create random reference sequence, and it's doubled version
    seq_len = 10000
    a_seq = utilities.create_r_seq(seq_len)
    d_seq = a_seq + a_seq

    # Create random candidate sequence
    r_seq = utilities.create_r_seq(seq_len)

    # Align random reference sequence to itself using all aligners
    a_seq_score_deflt_a = aligner_deflt.score(a_seq, a_seq)
    a_seq_score_jn_01_a = aligner_jn_01.score(a_seq, a_seq)
    a_seq_score_rl_01_a = aligner_rl_01.score(a_seq, a_seq)

    # Align random reference sequence to it's doubled version using
    # all aligners
    a_seq_score_deflt_d = aligner_deflt.score(a_seq, d_seq)
    a_seq_score_jn_01_d = aligner_jn_01.score(a_seq, d_seq)
    a_seq_score_rl_01_d = aligner_rl_01.score(a_seq, d_seq)

    # Align random candidate sequence to itself using all aligners
    r_seq_score_deflt_a = aligner_deflt.score(r_seq, a_seq)
    r_seq_score_jn_01_a = aligner_jn_01.score(r_seq, a_seq)
    r_seq_score_rl_01_a = aligner_rl_01.score(r_seq, a_seq)

    # Align random candidate sequence to it's doubled version using
    # all aligners
    r_seq_score_deflt_d = aligner_deflt.score(r_seq, d_seq)
    r_seq_score_jn_01_d = aligner_jn_01.score(r_seq, d_seq)
    r_seq_score_rl_01_d = aligner_rl_01.score(r_seq, d_seq)

    # Print summary results
    print("{0:8s} {1:8s} {2:8s} {3:8s} {4:8s}".format("", "actual", "", "random", ""))
    print(
        "{0:8s} {1:8s} {2:8s} {3:8s} {4:8s}".format(
            "aligner", "single", "double", "single", "double"
        )
    )
    print(
        "{0:8s} {1:8.1f} {2:8.1f} {3:8.1f} {4:8.1f}".format(
            "default",
            a_seq_score_deflt_a,
            a_seq_score_deflt_d,
            r_seq_score_deflt_a,
            r_seq_score_deflt_d,
        )
    )
    print(
        "{0:8s} {1:8.1f} {2:8.1f} {3:8.1f} {4:8.1f}".format(
            "jn_01",
            a_seq_score_jn_01_a,
            a_seq_score_jn_01_d,
            r_seq_score_jn_01_a,
            r_seq_score_jn_01_d,
        )
    )
    print(
        "{0:8s} {1:8.1f} {2:8.1f} {3:8.1f} {4:8.1f}".format(
            "rl_01",
            a_seq_score_rl_01_a,
            a_seq_score_rl_01_d,
            r_seq_score_rl_01_a,
            r_seq_score_rl_01_d,
        )
    )


def plot_alignment_with_random(reprocess=False):
    """Plot alignment of random candidate sequences of length 9000 to
    11,000 nt with a random reference sequence of length 10,000
    nt. The reference sequence is used as is, and doubled for
    alignment.

    Parameters
    ----------
    reprocess : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """
    # Read and parse configuration file
    config_dir = "../../resources"
    config_file = "pairwise-rl-01.cfg"
    config_rl_01 = ConfigParser()
    config_rl_01.read(os.path.join(config_dir, config_file))

    # Create aligner using the specified parameters
    aligner_rl_01 = utilities.create_aligner(config_rl_01["aligner"])

    # Process and dump, or load alignments
    pickle_file_name = config_file.replace(".cfg", "-random.pickle")
    if not os.path.exists(pickle_file_name) or reprocess:

        # Create random reference sequence, and it's doubled version
        seq_len = 10000
        a_seq = utilities.create_r_seq(seq_len)
        d_seq = a_seq + a_seq

        # Create and align random candidate sequences
        seq_lens = range(9000, 11000, 100)
        r_seq_score_a = []
        r_seq_score_d = []
        for seq_len in seq_lens:
            print("Aligning with length {:d}".format(seq_len))
            r_seq = utilities.create_r_seq(seq_len)
            r_seq_score_a.append(aligner_rl_01.score(r_seq, a_seq))
            r_seq_score_d.append(aligner_rl_01.score(r_seq, d_seq))

        seq_lens = np.array(seq_lens)
        r_seq_score_a = np.array(r_seq_score_a)
        r_seq_score_d = np.array(r_seq_score_d)

        # Dump alignments
        alignments = {}
        alignments["config"] = dict(config_rl_01)
        alignments["seq_len"] = seq_len
        alignments["seq_lens"] = seq_lens
        alignments["r_seq_score_a"] = r_seq_score_a
        alignments["r_seq_score_d"] = r_seq_score_d
        print("Dumping alignments")
        with open(pickle_file_name, "wb") as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, "rb") as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments["seq_len"]
        seq_lens = alignments["seq_lens"]
        r_seq_score_a = alignments["r_seq_score_a"]
        r_seq_score_d = alignments["r_seq_score_d"]

    # Plot and print absolute alignment score as a function of
    # candidate sequence length
    fig, ax = plt.subplots()
    ax.plot(seq_lens, r_seq_score_a, "r1")
    ax.plot(seq_lens, r_seq_score_d, "g2")
    ax.set_title("Reference sequence length: 10,000")
    ax.set_xlabel("Candidate sequence length")
    ax.set_ylabel("Absolute alignment score")
    plot_file_name = config_file.replace(".cfg", "-random.png")
    plt.savefig(plot_file_name)
    plt.show()

    return alignments


def plot_alignment_with_offsets(reprocess=False):
    """Plot alignment of random candidate and offset to reference
    sequences of length 10,000 nt. Offsets vary from 0 to 5,000 nt.
    The reference sequence is used as is, and doubled for alignment.

    Parameters
    ----------
    reprocess : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """
    # Read and parse configuration file
    config_dir = "../../resources"
    config_file = "pairwise-rl-01.cfg"
    config_rl_01 = ConfigParser()
    config_rl_01.read(os.path.join(config_dir, config_file))

    # Create aligner using the specified parameters
    aligner_rl_01 = utilities.create_aligner(config_rl_01["aligner"])

    # Process and dump, or load alignments
    pickle_file_name = config_file.replace(".cfg", "-offsets.pickle")
    if not os.path.exists(pickle_file_name) or reprocess:

        # Create random reference sequence, and it's doubled version
        seq_len = 10000
        a_seq = utilities.create_r_seq(seq_len)
        d_seq = a_seq + a_seq

        # Align offset to reference sequence, or its double
        o_seq_score_a = []
        o_seq_score_d = []
        n_offs = range(0, int(seq_len / 2), int(seq_len / 20))
        for n_off in n_offs:
            print("Aligning with offset {:d}".format(n_off))
            o_seq = create_o_seq(a_seq, n_off)
            o_seq_score_a.append(aligner_rl_01.score(o_seq, a_seq))
            o_seq_score_d.append(aligner_rl_01.score(o_seq, d_seq))

        n_offs = np.array(n_offs)
        o_seq_score_a = np.array(o_seq_score_a)
        o_seq_score_d = np.array(o_seq_score_d)

        # Dump alignments
        alignments = {}
        alignments["config"] = dict(config_rl_01)
        alignments["seq_len"] = seq_len
        alignments["n_offs"] = n_offs
        alignments["o_seq_score_a"] = o_seq_score_a
        alignments["o_seq_score_d"] = o_seq_score_d
        print("Dumping alignments")
        with open(pickle_file_name, "wb") as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, "rb") as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments["seq_len"]
        n_offs = alignments["n_offs"]
        o_seq_score_a = alignments["o_seq_score_a"]
        o_seq_score_d = alignments["o_seq_score_d"]

    # Plot and print absolute alignment score as a function of offset
    # length
    fig, ax = plt.subplots()
    ax.plot(n_offs, o_seq_score_a, "r")
    ax.plot(n_offs, o_seq_score_d, "g")
    ax.set_title("Reference sequence length 10,000")
    ax.set_xlabel("Number of nucleotides offset")
    ax.set_ylabel("Absolute alignment score")
    plot_file_name = config_file.replace(".cfg", "-offsets.png")
    plt.savefig(plot_file_name)
    plt.show()

    return alignments


def plot_alignment_with_errors(reprocess=False):
    """Plot alignment of error and reference sequences of length
    10,000 nt. Error sequences differ from the reference sequence by
    from 0 to 10,000 nts. The reference sequence is used as is, and
    doubled for alignment.

    Parameters
    ----------
    aligner_config : dct
        pairwise aligner configuration
    reprocess : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """
    # Read and parse configuration file
    config_dir = "../../resources"
    config_file = "pairwise-rl-01.cfg"
    config_rl_01 = ConfigParser()
    config_rl_01.read(os.path.join(config_dir, config_file))

    # Create aligner using the specified parameters
    aligner_rl_01 = utilities.create_aligner(config_rl_01["aligner"])

    # Process and dump, or load alignments
    pickle_file_name = config_file.replace(".cfg", "-errors.pickle")
    if not os.path.exists(pickle_file_name) or reprocess:

        # Create random reference sequence, and it's doubled version
        seq_len = 10000
        a_seq = utilities.create_r_seq(seq_len)
        d_seq = a_seq + a_seq

        # Create and align error candidate sequences
        e_seq_score_a = []
        e_seq_score_d = []
        n_errs = range(0, seq_len, int(seq_len / 50))
        for n_err in n_errs:
            print("Aligning with {:d} errors".format(n_err))
            e_seq = create_e_seq(a_seq, n_err)
            e_seq_score_a.append(aligner_rl_01.score(e_seq, a_seq))
            e_seq_score_d.append(aligner_rl_01.score(e_seq, d_seq))

        n_errs = np.array(n_errs)
        e_seq_score_a = np.array(e_seq_score_a)
        e_seq_score_d = np.array(e_seq_score_d)

        # Dump alignments
        alignments = {}
        alignments["config"] = dict(config_rl_01)
        alignments["seq_len"] = seq_len
        alignments["n_errs"] = n_errs
        alignments["e_seq_score_a"] = e_seq_score_a
        alignments["e_seq_score_d"] = e_seq_score_d
        print("Dumping alignments")
        with open(pickle_file_name, "wb") as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, "rb") as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments["seq_len"]
        n_errs = alignments["n_errs"]
        e_seq_score_a = alignments["e_seq_score_a"]
        e_seq_score_d = alignments["e_seq_score_d"]

    # Plot and print absolute alignment score as a function of error
    # length
    fig, ax = plt.subplots()
    ax.plot(n_errs, e_seq_score_a, "r--")
    ax.plot(n_errs, e_seq_score_d, "g:")
    ax.set_title("Reference sequence length 10,000")
    ax.set_xlabel("Number of nucleotides in error")
    ax.set_ylabel("Absolute alignment score")
    plot_file_name = config_file.replace(".cfg", "-errors.png")
    plt.savefig(plot_file_name)
    plt.show()

    return alignments


def plot_alignment_with_offsets_and_errors(reprocess=False):
    """Plot alignment of random candidate and offset with errors to
    reference sequences of length 10,000 nt. Offsets vary from 0 to
    5,000 nt.  Error sequences differ from the reference sequence by
    from 0 to 10,000 nts. The reference sequence is used as is, and
    doubled for alignment.

    Parameters
    ----------
    reprocess : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """
    # Read and parse configuration file
    config_dir = "../../resources"
    config_file = "pairwise-rl-01.cfg"
    config_rl_01 = ConfigParser()
    config_rl_01.read(os.path.join(config_dir, config_file))

    # Create aligner using the specified parameters
    aligner_rl_01 = utilities.create_aligner(config_rl_01["aligner"])

    # Process and dump, or load alignments
    pickle_file_name = config_file.replace(".cfg", "-offsets-and-errors.pickle")
    if not os.path.exists(pickle_file_name) or reprocess:

        # Create random reference sequence, and it's doubled version
        seq_len = 10000
        a_seq = utilities.create_r_seq(seq_len)
        d_seq = a_seq + a_seq

        # Create and align offset candidate sequences with errors
        n_offs = range(0, int(seq_len / 2), int(seq_len / 8))
        n_errs = range(0, int(seq_len / 2), int(seq_len / 20))

        o_e_seq_score_a = np.empty((len(n_offs), len(n_errs)))
        o_e_seq_score_d = np.empty((len(n_offs), len(n_errs)))
        r_o_e_seq_score_a = np.empty((len(n_offs), len(n_errs)))

        i_off = 0
        for n_off in n_offs:
            print("Aligning with offset {:d}".format(n_off))
            o_seq = create_o_seq(a_seq, n_off)

            i_err = 0
            for n_err in n_errs:
                print("Aligning with {:d} errors".format(n_err))
                o_e_seq = create_e_seq(o_seq, n_err)

                o_e_seq_score_a[i_off, i_err] = aligner_rl_01.score(o_e_seq, a_seq)
                o_e_seq_score_d[i_off, i_err] = aligner_rl_01.score(o_e_seq, d_seq)

                r_a_seq, r_o_e_seq = utilities.rotate_seqs(a_seq, o_e_seq)

                r_o_e_seq_score_a[i_off, i_err] = aligner_rl_01.score(
                    r_o_e_seq, r_a_seq
                )

                i_err += 1

            i_off += 1

        n_errs = np.array(n_errs)

        # Dump alignments
        alignments = {}
        alignments["config"] = dict(config_rl_01)
        alignments["seq_len"] = seq_len
        alignments["n_offs"] = n_offs
        alignments["n_errs"] = n_errs
        alignments["o_e_seq_score_a"] = o_e_seq_score_a
        alignments["o_e_seq_score_d"] = o_e_seq_score_d
        alignments["r_o_e_seq_score_a"] = r_o_e_seq_score_a
        print("Dumping alignments")
        with open(pickle_file_name, "wb") as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, "rb") as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments["seq_len"]
        n_offs = alignments["n_offs"]
        n_errs = alignments["n_errs"]
        o_e_seq_score_a = alignments["o_e_seq_score_a"]
        o_e_seq_score_d = alignments["o_e_seq_score_d"]
        r_o_e_seq_score_a = alignments["r_o_e_seq_score_a"]

    # Plot and print absolute alignment score as a function of error
    # length for each offset
    for i_fig in [1, 2, 3]:
        fig, ax = plt.subplots()

        i_off = 0
        for n_off in n_offs:
            ax.plot(n_errs, o_e_seq_score_a[i_off, :], "c")

            ax.plot([n_errs[0] + 100, 1900], o_e_seq_score_a[i_off, [0, 0]], "c:")
            ax.text(
                2000,
                o_e_seq_score_a[i_off, 0],
                n_off,
                verticalalignment="center",
                color="c",
            )

            if i_fig > 1:
                ax.plot(n_errs, r_o_e_seq_score_a[i_off, :], "b")

            if i_fig > 2:
                ax.plot(n_errs, o_e_seq_score_d[i_off, :], "k")

            i_off += 1

        ax.text(
            1000,
            o_e_seq_score_a[0, 0],
            "Offset",
            verticalalignment="center",
            color="c",
            bbox=dict(facecolor="w", edgecolor="w"),
        )

        if i_fig > 1:
            i_err = np.min(np.where(r_o_e_seq_score_a[0, :] < 5000))

            ax.plot(
                [n_errs[i_err] + 100, 2400], r_o_e_seq_score_a[0, [i_err, i_err]], "b:"
            )
            ax.text(
                2500,
                r_o_e_seq_score_a[0, i_err],
                "Optimal Rotation",
                verticalalignment="center",
                color="b",
            )

        if i_fig > 2:
            i_err = np.min(np.where(o_e_seq_score_d[0, :] < 3750))

            ax.plot(
                [n_errs[i_err] + 100, 2900], o_e_seq_score_a[0, [i_err, i_err]], "k:"
            )
            ax.text(
                3000,
                o_e_seq_score_d[0, i_err],
                "Doubled Reference",
                verticalalignment="center",
                color="k",
            )

        ax.set_title("Reference sequence length 10,000")
        ax.set_xlabel("Number of nucleotides in error")
        ax.set_ylabel("Absolute alignment score")

        plot_file_name = config_file.replace(
            ".cfg", "-offsets-and-errors-{}.png".format(i_fig)
        )
        plt.savefig(plot_file_name)
        plt.show()

    return alignments


if __name__ == "__main__":
    """Create some simple plots to understand operation of the
    Biopython PairwiseAligner.
    """
    # Add and parse arguments
    parser = ArgumentParser()

    parser.add_argument(
        "-n", "--plot-random", action="store_true", help="plot alignment with random"
    )
    parser.add_argument(
        "-o", "--plot-offsets", action="store_true", help="plot alignment with offsets"
    )
    parser.add_argument(
        "-e", "--plot-errors", action="store_true", help="plot alignment with errors"
    )
    parser.add_argument(
        "-m",
        "--plot-offsets-and-errors",
        action="store_true",
        help="plot alignment with offsets and errors",
    )

    parser.add_argument(
        "-r", "--reprocess", action="store_true", help="reprocess alignments"
    )

    parser.add_argument(
        "-a",
        "--test-aligners",
        action="store_true",
        help="test various aligner configurations",
    )
    parser.add_argument(
        "-t",
        "--test-rotation",
        action="store_true",
        help="test cyclic rotation of sequence",
    )

    options = parser.parse_args()

    # Plot alignment with random
    if options.plot_random:
        alignments_with_random = plot_alignment_with_random(options.reprocess)

    # Plot alignment with offsets
    if options.plot_offsets:
        alignments_with_offsets = plot_alignment_with_offsets(options.reprocess)

    # Plot alignment with errors
    if options.plot_errors:
        alignments_with_errors = plot_alignment_with_errors(options.reprocess)

    # Plot alignment with offsets and errors
    if options.plot_offsets_and_errors:
        alignments_with_offsets_and_errors = plot_alignment_with_offsets_and_errors(
            options.reprocess
        )

    # Test various aligner configurations
    if options.test_aligners:
        test_aligners()

    # Test cyclic rotation of sequence
    if options.test_rotation:

        # Create random reference sequence
        seq_len = 10000
        a_seq = utilities.create_r_seq(seq_len)

        # Offset random reference sequence
        n_off = int(seq_len / 2)
        o_seq = create_o_seq(a_seq, n_off)

        print(utilities.rotate_seqs(a_seq, o_seq))
