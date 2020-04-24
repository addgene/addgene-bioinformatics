from argparse import ArgumentParser
import os
import pickle
import random
import time

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np

import utilities


# TODO: Add error rate, outer distance standard deviation, indel
# fraction, indel extension probablility, and random seed:
# error_rate=0.020, standard_deviation=50, indel_fraction=0.15,
# indel_extended_prob=0.30, and random_seed=0, then move to utilities
def simulate_paired_reads_clean(seq, seq_nm, number_pairs=25000,
                                len_first_read=250, len_second_read=250,
                                outer_distance=500):
    rd1_seq_rcds = []
    rd2_seq_rcds = []
    seq_len = len(seq)
    for iP in range(number_pairs):

        i_rd1 = random.randint(0, seq_len - outer_distance)
        rd1_seq_rcd = SeqRecord(
            seq[i_rd1:i_rd1 + len_first_read],
            id=str(iP), name="one", description="first read of read pair")
        rd1_seq_rcd.letter_annotations["phred_quality"] = (
            40 * np.ones(len_first_read, dtype=int)).tolist()
        rd1_seq_rcds.append(rd1_seq_rcd)

        i_rd2 = i_rd1 + outer_distance - len_second_read
        rd2_seq_rcd = SeqRecord(
            seq[i_rd2:i_rd2 + len_first_read],
            id=str(iP), name="two", description="second read of read pair")
        rd2_seq_rcd.letter_annotations["phred_quality"] = (
            40 * np.ones(len_second_read, dtype=int)).tolist()
        rd2_seq_rcds.append(rd2_seq_rcd)

    rd1_fNm = seq_nm + "_r1.fastq"
    rd2_fNm = seq_nm + "_r2.fastq"
    SeqIO.write(rd1_seq_rcds, rd1_fNm, "fastq")
    SeqIO.write(rd2_seq_rcds, rd2_fNm, "fastq")
    return rd1_fNm, rd2_fNm


# TODO: Complete in order to validate use of wgsim
def simulate_paired_reads_wgsim(seq, seq_nm, number_pairs=25000,
                                len_first_read=250, len_second_read=250,
                                outer_distance=500):
    seq_rcd = SeqRecord(seq, id="0", name="base", description="reference")
    seq_fNm = seq_nm + ".fasta"
    rd1_fNm = seq_nm + "_r1.fastq"
    rd2_fNm = seq_nm + "_r2.fastq"
    SeqIO.write(seq_rcd, seq_fNm, "fasta")
    utilities.wgsim(seq_fNm,
                    rd1_fNm,
                    rd2_fNm,
                    error_rate=0.0,
                    outer_distance=outer_distance,
                    standard_deviation=0,
                    number_pairs=number_pairs,
                    len_first_read=len_first_read,
                    len_second_read=len_second_read,
                    mutation_rate=0.0,
                    indel_fraction=0.0,
                    indel_extended_prob=0.0,
                    random_seed=0,
                    ambiguous_base_frac=0.0,
                    haplotype_mode=True)
    return rd1_fNm, rd2_fNm


if __name__ == "__main__":
    # Add and parse arguments
    parser = ArgumentParser()

    parser.add_argument('-w', '--with-repeats',
                        action='store_true',
                        help="create a random sequence with repeats")

    parser.add_argument('-p', '--plot-histogram',
                        action='store_true',
                        help="plot histogram of k-mer counts")

    parser.add_argument('-r', '--reprocess',
                        action='store_true',
                        help="reprocess alignments")

    options = parser.parse_args()

    pickle_file_name = __name__+"pickle"
    if not os.path.exists(pickle_file_name) or options.reprocess:
        start_time = time.time()
        print("creating random sequence with repeats ...", end=" ", flush=True)
        if options.with_repeats:
            r_seq = utilities.create_r_seq_w_rep()
        else:
            r_seq = utilities.create_r_seq(10000)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("random sequence length: {0}".format(len(r_seq)))

        start_time = time.time()
        print("counting k_mers in sequence ...", end=" ", flush=True)
        k_mers_in_seq = utilities.count_k_mers_in_seq(r_seq + r_seq)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # === clean ===

        start_time = time.time()
        print("simulating paired reads ...", end=" ", flush=True)
        rd1_fNm, rd2_fNm = simulate_paired_reads_clean(
            r_seq + r_seq, "ccc_clean")
        print("done in {0} s".format(time.time() - start_time), flush=True)

        start_time = time.time()
        print("counting k_mers in reads ...", end=" ", flush=True)
        k_mers_in_rd1_c, seq_rcds_rd1_c = utilities.count_k_mers_in_rds(rd1_fNm)
        k_mers_in_rd2_c, seq_rcds_rd2_c = utilities.count_k_mers_in_rds(rd2_fNm)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        start_time = time.time()
        print("writing k_mers and counts in reads ...", end=" ", flush=True)
        utilities.write_k_mer_counts_in_rds(
            k_mers_in_rd1_c, k_mers_in_rd2_c, "ccc_clean_cnt.txt")
        print("done in {0} s".format(time.time() - start_time), flush=True)

        start_time = time.time()
        print("counting k_mers in reads using kmc ...", end=" ", flush=True)
        inp_read_file_names = [rd1_fNm, rd2_fNm]
        inp_database_file_name = "ccc_clean"
        out_database_file_name = "ccc_clean_kmc.txt"
        utilities.kmc(inp_read_file_names, inp_database_file_name,
                      count_min=0, max_count=1e9, canonical_form=False)
        utilities.kmc_transform(inp_database_file_name, "dump",
                                out_database_file_name)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        start_time = time.time()
        print("reading k_mers counted using kmc ...", end=" ", flush=True)
        k_mers_by_kmc_c = utilities.read_k_mer_counts(out_database_file_name)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        cn_in_seq = np.array(
            [val['cnt'] / 2.0 for val in k_mers_in_seq.values()]
        )
        cn_in_rds_c = np.array(
            [(k_mers_in_rd1_c[key]['cnt'] + k_mers_in_rd2_c[key]['cnt']) / 2.0
             for key in k_mers_in_rd1_c.keys()]
        )
        cn_by_kmc_c = np.array(
            [val['cnt'] / 2.0 for val in k_mers_by_kmc_c.values()]
        )

        coverage = int(np.sum(cn_by_kmc_c) / np.sum(cn_in_seq))

        counts = {}
        counts['r_seq'] = r_seq
        counts['k_mers_in_seq'] = k_mers_in_seq
        counts['k_mers_in_rd1_c'] = k_mers_in_rd1_c
        counts['k_mers_in_rd2_c'] = k_mers_in_rd2_c
        counts['k_mers_by_kmc_c'] = k_mers_by_kmc_c
        counts['cn_in_seq'] = cn_in_seq
        counts['cn_in_rds_c'] = cn_in_rds_c
        counts['cn_by_kmc_c'] = cn_by_kmc_c
        counts['coverage'] = coverage

        print("Dumping counts")
        with open(pickle_file_name, 'wb') as pickle_file:
            pickle.dump(counts, pickle_file)

    else:

        print("Loading counts")
        with open(pickle_file_name, 'rb') as pickle_file:
            counts = pickle.load(pickle_file)

        r_seq = counts['r_seq']
        k_mers_in_seq = counts['k_mers_in_seq']
        k_mers_in_rd1_c = counts['k_mers_in_rd1_c']
        k_mers_in_rd2_c = counts['k_mers_in_rd2_c']
        k_mers_by_kmc_c = counts['k_mers_by_kmc_c']
        cn_in_seq = counts['cn_in_seq']
        cn_in_rds_c = counts['cn_in_rds_c']
        cn_by_kmc_c = counts['cn_by_kmc_c']
        coverage = counts['coverage']

    print("r_seq: ", r_seq)
    print("r_seq + r_seq: ", r_seq + r_seq)
    for k_mer in k_mers_in_seq.keys():
        ln = "k_mer: {0} in_seq: {1:6d}".format(
            k_mer, k_mers_in_seq[k_mer]["cnt"])
        if k_mer in k_mers_in_rd1_c.keys():
            ln += ", in_rd1_c: {0:6d}".format(k_mers_in_rd1_c[k_mer]["cnt"])
        else:
            ln += ", in_rd1_c: {0:6d}".format(0)
        if k_mer in k_mers_in_rd2_c.keys():
            ln += ", in_rd2_c: {0:6d}".format(k_mers_in_rd2_c[k_mer]["cnt"])
        else:
            ln += ", in_rd2_c: {0:6d}".format(0)
        if k_mer in k_mers_by_kmc_c.keys():
            ln += ", by_kmc_c: {0:6d}".format(k_mers_by_kmc_c[k_mer]["cnt"])
        else:
            ln += ", by_kmc_c: {0:6d}".format(0)
        print(ln)

    if options.plot_histogram:

        cn_in_rds_c = cn_in_rds_c / coverage
        cn_by_kmc_c = cn_by_kmc_c / coverage

        bin_stp = 1.0
        bin_min = min(
            np.min(cn_in_seq), np.min(cn_in_rds_c), np.min(cn_by_kmc_c)
        ) - bin_stp / 2
        bin_max = max(
            np.max(cn_in_seq), np.max(cn_in_rds_c), np.max(cn_by_kmc_c)
        ) + bin_stp / 2

        bin_edges = np.arange(bin_min, bin_max, bin_stp)

        cn_in_seq_hist, _ = np.histogram(cn_in_seq, bin_edges)

        cn_in_rds_c_hist, _ = np.histogram(cn_in_rds_c, bin_edges)
        cn_by_kmc_c_hist, _ = np.histogram(cn_by_kmc_c, bin_edges)

        fig, ax = plt.subplots()

        ax.bar(bin_edges[:-1] + bin_stp / 2 + 1 * bin_stp / 6, cn_in_seq_hist,
               width=0.3, align='center', color='r')

        ax.bar(bin_edges[:-1] + bin_stp / 2 + 3 * bin_stp / 6,
               cn_in_rds_c_hist,
               width=0.3, align='center', color='g')

        ax.bar(bin_edges[:-1] + bin_stp / 2 + 5 * bin_stp / 6,
               cn_by_kmc_c_hist,
               width=0.3, align='center', color='b')

        ax.set_title("Random sequence with repeats")
        ax.set_xlabel("k-mer counts")
        ax.set_ylabel("Frequency")

        x1, x2 = ax.get_xlim()
        xW = x2 - x1

        y1, y2 = ax.get_ylim()
        yW = y2 - y1

        ax.text(x1 + 0.10 * xW, y2 - 0.05 * yW,
                "Coverage: {}".format(coverage))

        plt.show()
