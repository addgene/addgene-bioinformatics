from argparse import ArgumentParser
import os
import pickle
from random import choice, choices

from utilities import create_r_seq

from Bio import Align
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import matplotlib.pyplot as plt
import numpy as np


def create_o_seq(a_seq, n_off):
    """Offsets a sequence by the specified number of nucleotides.

    Parameters
    ----------
    a_seq : Bio.Seq.Seq
        the seqeunce to offset
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
    """
    seq_len = len(o_seq)
    e_seq_lst = list(o_seq)
    for nt_idx in choices(range(seq_len), k=n_err):
        nt_lst = list("CGTA")
        del(nt_lst[nt_lst.index(e_seq_lst[nt_idx])])
        e_seq_lst[nt_idx] = choice(nt_lst)
    e_seq_str = "".join(e_seq_lst)
    e_seq = Seq(e_seq_str, IUPAC.unambiguous_dna)
    return e_seq


def align(e_seq, a_seq, mode='global',
          match_score=1.0, mismatch_score=0.0, gap_score=0.0):
    """Aligns the specified sequence using the specified parameters.

    Parameters
    ----------
    e_seq : Bio.Seq.Seq
        the sequence thought to contain errors, although any sequence
        is valid
    a_seq : Bio.Seq.Seq
        the sequence thought to represent truth, although any sequence
        is valid
    mode " str
        the alignment mode: 'global' | 'local'
    match_score : flt
        the score of a match
    mismatch_score : flt
        the score of a mismatch
    gap_score : flt
        the score of a gap

    Returns
    -------
    flt
        the score of the alignment
    """
    # Note BLAST defaults:
    #     aligner.match_score = 1.0
    #     aligner.mismatch_score = -2.0
    #     aligner.gap_score = -2.5
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.gap_score = gap_score
    score = aligner.score(a_seq, e_seq)
    return score


def plot_alignment_with_random(force=False):
    """Plot alignment of random candidate sequences of length 1000 to
    30,000 nt with a random reference sequence of length 10,000
    nt. The reference sequeence is doubled for alignment.

    Parameters
    ----------
    force : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """

    # Process and dump, or load alignments
    pickle_file_name = "plot_alignment_with_random.pickle"
    if not os.path.exists(pickle_file_name) or force:

        # Create random reference sequence
        seq_len = 10000
        c_seq = create_r_seq(seq_len)
        a_seq = c_seq + c_seq

        # Create and align random candidate sequences
        seq_lens = range(1000, 30000, 500)
        r_seq_score = []
        for seq_len in seq_lens:
            print("Aligning with length {:d}".format(seq_len))
            r_seq = create_r_seq(seq_len)
            r_seq_score.append(align(r_seq, a_seq))

        seq_lens = np.array(seq_lens)
        r_seq_score = np.array(r_seq_score)

        # Dump alignments
        alignments = {}
        alignments['seq_len'] = seq_len
        alignments['seq_lens'] = seq_lens
        alignments['r_seq_score'] = r_seq_score
        print("Dumping alignments")
        with open(pickle_file_name, 'wb') as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, 'rb') as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments['seq_len']
        seq_lens = alignments['seq_lens']
        r_seq_score = alignments['r_seq_score']

    # Plot absolute alignment score as a function of candidate
    # sequence length
    fig, ax = plt.subplots()
    ax.plot(seq_lens, r_seq_score, 'r-')
    ax.set_title("Reference sequence length: 10,000")
    ax.set_xlabel("Candidate sequence length")
    ax.set_ylabel("Absolute alignment score")
    plt.show()

    # Plot relative score as a function of candidate sequence length
    fig, ax = plt.subplots()
    ax.plot(seq_lens, r_seq_score / seq_lens, 'r-')
    ax.set_title("Reference sequence length: 10,000")
    ax.set_xlabel("Candidate sequence length")
    ax.set_ylabel("Relative alignment score")
    plt.show()

    return alignments


def plot_alignment_with_offsets(force=False):
    """Plot alignment of random candidate and offset to reference
    sequences of length 10,000 nt. Offsets vary from 0 to 5,000 nt.
    The reference sequeence is doubled for alignment, and global and
    local alignment modes are used.

    Parameters
    ----------
    force : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """

    # Process and dump, or load alignments
    pickle_file_name = "plot_alignment_with_offsets.pickle"
    if not os.path.exists(pickle_file_name) or force:

        # Create random reference sequence
        seq_len = 10000
        a_seq_1 = create_r_seq(seq_len)
        a_seq_2 = a_seq_1 + a_seq_1

        # Create random candidate sequence
        r_seq = create_r_seq(seq_len)

        # Align candidate to reference sequence, or its double, using
        # global and local alignment
        r_seq_score_1_g = np.array(align(r_seq, a_seq_1, mode='global'))
        r_seq_score_1_l = np.array(align(r_seq, a_seq_1, mode='local'))
        r_seq_score_2_g = np.array(align(r_seq, a_seq_2, mode='global'))
        r_seq_score_2_l = np.array(align(r_seq, a_seq_2, mode='local'))

        # Align offset to reference sequence, or its double, using
        # global and local alignment
        o_seq_score_1_g = []
        o_seq_score_1_l = []
        o_seq_score_2_g = []
        o_seq_score_2_l = []
        n_offs = range(0, int(seq_len / 2), int(seq_len / 20))
        for n_off in n_offs:
            print("Aligning with offset {:d}".format(n_off))
            o_seq = create_o_seq(a_seq_1, n_off)
            o_seq_score_1_g.append(align(o_seq, a_seq_1, mode='global'))
            o_seq_score_1_l.append(align(o_seq, a_seq_1, mode='local'))
            o_seq_score_2_g.append(align(o_seq, a_seq_2, mode='global'))
            o_seq_score_2_l.append(align(o_seq, a_seq_2, mode='local'))

        n_offs = np.array(n_offs)
        o_seq_score_1_g = np.array(o_seq_score_1_g)
        o_seq_score_1_l = np.array(o_seq_score_1_l)
        o_seq_score_2_g = np.array(o_seq_score_2_g)
        o_seq_score_2_l = np.array(o_seq_score_2_l)

        # Dump alignments
        alignments = {}
        alignments['seq_len'] = seq_len
        alignments['r_seq_score_1_g'] = r_seq_score_1_g
        alignments['r_seq_score_1_l'] = r_seq_score_1_l
        alignments['r_seq_score_2_g'] = r_seq_score_2_g
        alignments['r_seq_score_2_l'] = r_seq_score_2_l
        alignments['n_offs'] = n_offs
        alignments['o_seq_score_1_g'] = o_seq_score_1_g
        alignments['o_seq_score_1_l'] = o_seq_score_1_l
        alignments['o_seq_score_2_g'] = o_seq_score_2_g
        alignments['o_seq_score_2_l'] = o_seq_score_2_l
        print("Dumping alignments")
        with open(pickle_file_name, 'wb') as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, 'rb') as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments['seq_len']
        r_seq_score_1_g = alignments['r_seq_score_1_g']
        r_seq_score_1_l = alignments['r_seq_score_1_l']
        r_seq_score_2_g = alignments['r_seq_score_2_g']
        r_seq_score_2_l = alignments['r_seq_score_2_l']
        n_offs = alignments['n_offs']
        o_seq_score_1_g = alignments['o_seq_score_1_g']
        o_seq_score_1_l = alignments['o_seq_score_1_l']
        o_seq_score_2_g = alignments['o_seq_score_2_g']
        o_seq_score_2_l = alignments['o_seq_score_2_l']

    # Plot absolute alignment score as a function of offset length
    fig, ax = plt.subplots()
    ax.plot([min(n_offs), max(n_offs)],
            [r_seq_score_1_g, r_seq_score_1_g], 'b-')
    ax.plot([min(n_offs), max(n_offs)],
            [r_seq_score_1_l, r_seq_score_1_l], 'bo')
    ax.plot([min(n_offs), max(n_offs)],
            [r_seq_score_2_g, r_seq_score_2_g], 'k-')
    ax.plot([min(n_offs), max(n_offs)],
            [r_seq_score_2_l, r_seq_score_2_l], 'ko')
    ax.plot(n_offs, o_seq_score_1_g, 'r-')
    ax.plot(n_offs, o_seq_score_1_l, 'ro')
    ax.plot(n_offs, o_seq_score_2_g, 'g-')
    ax.plot(n_offs, o_seq_score_2_l, 'go')
    ax.set_title("Reference sequence length 10,000")
    ax.set_xlabel("Number of nucleotides offset")
    ax.set_ylabel("Absolute alignment score")
    plt.show()

    return alignments


def plot_alignment_with_errors(force=False):
    """Plot alignment of error and reference sequences of length
    10,000 nt. Error sequences differ from the reference sequence by
    from 0 to 10,000 nts.  The reference sequeence is doubled for
    alignment.

    Parameters
    ----------
    force : bln
        flag to reprocess alignments, or not

    Returns
    -------
    dct
        dictionary of alignment values
    """

    # Process and dump, or load alignments
    pickle_file_name = "plot_alignment_with_errors.pickle"
    if not os.path.exists(pickle_file_name) or force:

        # Create random reference sequence
        seq_len = 10000
        c_seq = create_r_seq(seq_len)
        a_seq = c_seq + c_seq

        # Create and align random candidate sequence
        r_seq = create_r_seq(seq_len)
        r_seq_score = align(r_seq, a_seq)

        # Create and align error candidate sequences
        e_seq_score = []
        n_errs = range(0, seq_len, int(seq_len / 50))
        for n_err in n_errs:
            print("Aligning with {:d} errors".format(n_err))
            e_seq = create_e_seq(c_seq, n_err)
            e_seq_score.append(align(e_seq, a_seq))

        r_seq_score = np.array(r_seq_score)
        n_errs = np.array(n_errs)
        e_seq_score = np.array(e_seq_score)

        # Dump alignments
        alignments = {}
        alignments['seq_len'] = seq_len
        alignments['r_seq_score'] = r_seq_score
        alignments['n_errs'] = n_errs
        alignments['e_seq_score'] = e_seq_score
        print("Dumping alignments")
        with open(pickle_file_name, 'wb') as pickle_file:
            pickle.dump(alignments, pickle_file)

    else:

        # Load alignments
        print("Loading alignments")
        with open(pickle_file_name, 'rb') as pickle_file:
            alignments = pickle.load(pickle_file)
        seq_len = alignments['seq_len']
        r_seq_score = alignments['r_seq_score']
        n_errs = alignments['n_errs']
        e_seq_score = alignments['e_seq_score']

    # Plot absolute alignment score as a function of error length
    fig, ax = plt.subplots()
    ax.plot([min(n_errs), max(n_errs)],
            [r_seq_score, r_seq_score], 'b-')
    ax.plot(n_errs, e_seq_score, 'r-')
    ax.set_title("Reference sequence length 10,000")
    ax.set_xlabel("Number of nucleotides in error")
    ax.set_ylabel("Absolute alignment score")
    plt.show()

    # Plot relative alignment score as a function of error length
    fig, ax = plt.subplots()
    rel_n_errs = 100.0 * n_errs / seq_len
    rel_e_seq_score = (
        100.0 * (e_seq_score - r_seq_score) / (seq_len - r_seq_score))
    ax.plot(rel_n_errs, rel_e_seq_score, 'r-')
    ax.set_title("Reference sequence length 10,000")
    ax.set_xlabel("Number of nucleotides in error")
    ax.set_ylabel("Relative alignment score")
    plt.show()

    return alignments


if __name__ == "__main__":
    """Create some simple plots to understand operation of the
    Biopython PairwiseAligner.
    """

    # And and parse arguments
    parser = ArgumentParser()
    parser.add_argument('-m', '--plot-random',
                        action='store_true',
                        help="plot alignment with random")
    parser.add_argument('-o', '--plot-offsets',
                        action='store_true',
                        help="plot alignment with offsets")
    parser.add_argument('-e', '--plot-errors',
                        action='store_true',
                        help="plot alignment with errors")
    parser.add_argument('-r', '--reprocess',
                        action='store_true',
                        help="reprocess alignments")
    options = parser.parse_args()

    # Plot alignment with random
    if options.plot_random:
        alignments_with_random = plot_alignment_with_random(
            options.reprocess)

    # Plot alignment with offsets
    if options.plot_offsets:
        alignments_with_offsets = plot_alignment_with_offsets(
            options.reprocess)

    # Plot alignment with errors
    if options.plot_errors:
        alignments_with_errors = plot_alignment_with_errors(
            options.reprocess)
