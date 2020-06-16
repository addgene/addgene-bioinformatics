from argparse import ArgumentParser
from configparser import ConfigParser
import os
import pickle
import shutil
import time

from Bio import SeqIO
import numpy as np
from sklearn.cluster import KMeans

import utilities


# Set processing parameters
BASE_FILE_NAME = "random_seq"
SEQ_LEN = 10000
K_MER_LEN = 25
EXP_CNT = 16
K_MER_CNT_REP = [EXP_CNT]
NUMBER_PAIRS = 25000
LEN_FIRST_READ = 250
LEN_SECOND_READ = 250
OUTER_DISTANCE = 500
FRAGMENT_LEN = OUTER_DISTANCE
LEN_UNPAIRED_READ = 500
NUMBER_READS = 250


def collect_k_mer_cnt(k_mers, seq_cnt=2):
    """Collect k-mer counts, assuming that we counted in a doubled
    sequence. Set sequence count otherwise.
    """
    return np.array(
        [val['cnt'] / seq_cnt for val in k_mers.values()]
    )


def write_paired_reads_for_cnt(
        k_mer_cnt_rd1, coverage_rd1, k_mers_rd1, seq_rcds_rd1,
        k_mer_cnt_rd2, coverage_rd2, k_mers_rd2, seq_rcds_rd2,
        K_MER_CNT_REP, EXP_CNT):

    if len(seq_rcds_rd1) != len(seq_rcds_rd2):
        raise(Exception("Number of sequence records for each pair must agree"))

    i_rcds_rd1, _ = find_seq_rcds_for_cnt(
        k_mer_cnt_rd1, coverage_rd1, k_mers_rd1)
    i_rcds_rd2, _ = find_seq_rcds_for_cnt(
        k_mer_cnt_rd2, coverage_rd2, k_mers_rd2)
    i_rcds = i_rcds_rd1.intersection(i_rcds_rd2)

    rd1_wr_file_name, rd1_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd1, "rd1")
    rd2_wr_file_name, rd2_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd2, "rd2")

    return (rd1_wr_file_name, rd1_wo_file_name,
            rd2_wr_file_name, rd2_wo_file_name)


def find_seq_rcds_for_cnt(k_mer_cnt_rds, coverage, k_mers):
    """Use k-means clustering to identify reads containing a k-mer
    with the expected count.
    """
    # Cluster counts in the expected number of clusters
    kmeans = KMeans(n_clusters=len(K_MER_CNT_REP) + 1, random_state=0).fit(
        (k_mer_cnt_rds / coverage).reshape(-1, 1)
    )

    # Identify the cluster and sequence records (reads) corresponding
    # to the expected count
    label = kmeans.predict(np.array([EXP_CNT]).reshape(-1, 1))[0]
    i_rcds = set()
    keys = np.array(list(k_mers.keys()))
    for key in keys[kmeans.labels_ == label]:
        i_rcds = i_rcds.union(k_mers[key]['src'])

    return i_rcds, kmeans


def write_reads_for_cnt(i_rcds, seq_rcds, case):
    """Write a FASTQ file containing the sequence records (reads)
    specified, and a file containing all other sequence records
    (reads).
    """
    # Write a FASTQ file containing the specified sequence records
    # (reads)
    read_wr_file_name = BASE_FILE_NAME + "_wr_" + case + ".fastq"
    with open(read_wr_file_name, "w") as f:
        for i_rcd in i_rcds:
            SeqIO.write(seq_rcds[i_rcd], f, "fastq")

    # Write a FASTQ file containing all other sequence records (reads)
    read_wo_file_name = BASE_FILE_NAME + "_wo_" + case + ".fastq"
    with open(read_wo_file_name, "w") as f:
        for i_rcd in set(range(NUMBER_PAIRS)).difference(i_rcds):
            SeqIO.write(seq_rcds[i_rcd], f, "fastq")

    return read_wr_file_name, read_wo_file_name


if __name__ == "__main__":

    # Add and parser arguments
    parser = ArgumentParser()
    parser.add_argument('-i', '--initialize',
                        action='store_true',
                        help="create sequence and count k-mers")
    parser.add_argument('-r', '--with-repeats',
                        action='store_true',
                        help="create and process a random sequence with repeats")
    parser.add_argument('-s', '--assemble-using-spades',
                        action='store_true',
                        help="assemble paired reads using SPAdes")
    parser.add_argument('-u', '--assemble-using-unicycler',
                        action='store_true',
                        help="assemble paired reads using Unicycler")
    options = parser.parse_args()

    # Read and parse aligner configuration file
    home_dir = os.path.join(
        os.sep, *os.path.abspath(__file__).split(os.sep)[0:-4])
    aligner_config_dir = os.path.join(home_dir, "resources")
    aligner_config_file = "pairwise-rl-01.cfg"
    aligner_config = ConfigParser()
    aligner_config.read(os.path.join(aligner_config_dir, aligner_config_file))
    aligner_config['aligner']['file'] = aligner_config_file

    # Create a pairwise aligner with the specified configuration
    aligner = utilities.create_aligner(aligner_config['aligner'])

    # Create and dump, or load results
    if options.with_repeats:
        BASE_FILE_NAME = "repetitive_seq"
    pickle_file_name = BASE_FILE_NAME + ".pickle"
    if not os.path.exists(pickle_file_name) or options.initialize:

        # Create a random sequence with repeats, or if not, a random
        # sequence of the same length
        start_time = time.time()
        if options.with_repeats:
            print("Creating random sequence with repeats ...", end=" ",
                  flush=True)
            # Note that repeats need to occur in the interior of the
            # sequence in order to obtain the expected count since we
            # are simapling in a linear sequence. This issue may be
            # resolved with the random rotation added to the paired
            # read simulation
            a_seq = utilities.create_r_seq(
                OUTER_DISTANCE
            )
            b_seq, r_k_mers = utilities.create_r_seq_w_rep(
                k_mer_len=K_MER_LEN, k_mer_cnt=K_MER_CNT_REP
            )
            c_seq = utilities.create_r_seq(
                SEQ_LEN - OUTER_DISTANCE - sum(K_MER_CNT_REP) * K_MER_LEN
            )
            r_seq = a_seq + b_seq + c_seq
        else:
            print("Creating random sequence without repeats ...", end=" ",
                  flush=True)
            r_seq = utilities.create_r_seq(SEQ_LEN)
            r_k_mers = []
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Random sequence length: {0}".format(len(r_seq)))

        # Count k-mers in the random sequence, doubled to represent
        # circular DNA
        start_time = time.time()
        print("Counting k_mers in sequence ...", end=" ", flush=True)
        k_mers_seq = utilities.count_k_mers_in_seq(r_seq + r_seq)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # Simulate paired reads of the random sequence
        start_time = time.time()
        print("Simulating paired reads from random sequence ...", end=" ",
              flush=True)
        rd1_file_name, rd2_file_name = utilities.simulate_paired_reads(
            r_seq,
            BASE_FILE_NAME,
            number_pairs=NUMBER_PAIRS,
            len_first_read=LEN_FIRST_READ,
            len_second_read=LEN_SECOND_READ,
            outer_distance=OUTER_DISTANCE)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # Count k-mers in the paired reads, and write the result to a
        # file
        start_time = time.time()
        print("Counting k_mers in reads ...", end=" ", flush=True)
        k_mers_rd1, seq_rcds_rd1 = utilities.count_k_mers_in_rds(
            rd1_file_name)
        k_mers_rd2, seq_rcds_rd2 = utilities.count_k_mers_in_rds(
            rd2_file_name)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        start_time = time.time()
        print("Writing k_mers and counts in reads ...", end=" ", flush=True)
        utilities.write_k_mer_counts_in_rds(
            k_mers_rd1, k_mers_rd2, BASE_FILE_NAME + "_cnt.txt")
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # Create the results dictionary, and dump to a pickle file
        start_time = time.time()
        print("Dumping results ...", end=" ", flush=True)
        results = {}
        results['r_seq'] = r_seq
        results['r_k_mers'] = r_k_mers
        results['rd1_file_name'] = rd1_file_name
        results['rd2_file_name'] = rd2_file_name
        results['k_mers_seq'] = k_mers_seq
        results['k_mers_rd1'] = k_mers_rd1
        results['k_mers_rd2'] = k_mers_rd2
        results['seq_rcds_rd1'] = seq_rcds_rd1
        results['seq_rcds_rd2'] = seq_rcds_rd2
        with open(pickle_file_name, 'wb') as pickle_file:
            pickle.dump(results, pickle_file)
        print("done in {0} s".format(time.time() - start_time), flush=True)

    else:

        # Load a pickle file, and extract values from the results
        # dictionary
        start_time = time.time()
        print("Loading results ...", end=" ", flush=True)
        with open(pickle_file_name, 'rb') as pickle_file:
            results = pickle.load(pickle_file)
        r_seq = results['r_seq']
        r_k_mers = results['r_k_mers']
        rd1_file_name = results['rd1_file_name']
        rd2_file_name = results['rd2_file_name']
        k_mers_seq = results['k_mers_seq']
        k_mers_rd1 = results['k_mers_rd1']
        k_mers_rd2 = results['k_mers_rd2']
        seq_rcds_rd1 = results['seq_rcds_rd1']
        seq_rcds_rd2 = results['seq_rcds_rd2']
        print("done in {0} s".format(time.time() - start_time), flush=True)

    # Perform processing for random sequences
    if not options.with_repeats:

        # Assemble paired reads using SPAdes
        spades_wo_out_dir = BASE_FILE_NAME + "_spades_wo"
        if options.assemble_using_spades:
            start_time = time.time()
            print("Assembling paired reads using SPAdes ...", end=" ",
                  flush=True)
            if os.path.exists(spades_wo_out_dir):
                shutil.rmtree(spades_wo_out_dir)
            command = utilities.spades(rd1_file_name,
                                       rd2_file_name,
                                       spades_wo_out_dir)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)

        # Read the multi-FASTA SPAdes output file, and align to the
        # doubled random sequence
        if os.path.exists(spades_wo_out_dir):
            spades_seq_rcds = [seq_record for seq_record in SeqIO.parse(
                os.path.join(spades_wo_out_dir,
                             "contigs.fasta"), "fasta")]
            print("SPAdes score: {0}".format(
                aligner.score(r_seq + r_seq, spades_seq_rcds[0].seq)))

        # Assemble paired reads using Unicycler
        unicycler_wo_out_dir = BASE_FILE_NAME + "_unicycler_wo"
        if options.assemble_using_unicycler:
            start_time = time.time()
            print("Assembling paired reads using Unicycler ...", end=" ",
                  flush=True)
            if os.path.exists(unicycler_wo_out_dir):
                shutil.rmtree(unicycler_wo_out_dir)
            command = utilities.unicycler(rd1_file_name,
                                          rd2_file_name,
                                          unicycler_wo_out_dir)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)

        # Read the multi-FASTA Unicycler output file, and align to the
        # doubled random sequence
        if os.path.exists(unicycler_wo_out_dir):
            unicycler_seq_rcds = [seq_record for seq_record in SeqIO.parse(
                os.path.join(unicycler_wo_out_dir,
                             "assembly.fasta"), "fasta")]
            print("Unicycler score: {0}".format(
                aligner.score(r_seq + r_seq, unicycler_seq_rcds[0].seq)))

        pass

    # Perform processing for repetitive sequences
    else:

        # Collect k-mer counts, and compute coverage
        k_mer_cnt_seq = collect_k_mer_cnt(k_mers_seq, seq_cnt=2)
        k_mer_cnt_rd1 = collect_k_mer_cnt(k_mers_rd1, seq_cnt=1)
        k_mer_cnt_rd2 = collect_k_mer_cnt(k_mers_rd2, seq_cnt=1)
        coverage_rd1 = int(np.sum(k_mer_cnt_rd1) / np.sum(k_mer_cnt_seq))
        coverage_rd2 = int(np.sum(k_mer_cnt_rd2) / np.sum(k_mer_cnt_seq))

        # Separate reads into those containing a k-mer with the
        # expected count, and all others
        print("Writing reads based on read count ...", end=" ", flush=True)
        r_k_mer_exp_rd1 = str(r_k_mers[K_MER_CNT_REP.index(EXP_CNT)])
        r_k_mer_exp_rd2 = r_k_mer_exp_rd1[-1::-1]
        (rd1_wr_file_name, rd1_wo_file_name,
         rd2_wr_file_name, rd2_wo_file_name) = write_paired_reads_for_cnt(
             k_mer_cnt_rd1, coverage_rd1, k_mers_rd1, seq_rcds_rd1,
             k_mer_cnt_rd2, coverage_rd2, k_mers_rd2, seq_rcds_rd2,
             K_MER_CNT_REP, EXP_CNT)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("r_k_mer_exp_rd1: {0}".format(r_k_mer_exp_rd1))
        print("r_k_mer_exp_rd2: {0}".format(r_k_mer_exp_rd2))

        # Assemble paired reads with repeats using SSAKE
        start_time = time.time()
        print("Assembling paired reads with repeats using SSAKE ...",
              end=" ", flush=True)
        ssake_out_base_name = BASE_FILE_NAME + "_ssake"
        if os.path.exists(ssake_out_base_name):
            shutil.rmtree(ssake_out_base_name)
        utilities.ssake(rd1_wr_file_name,
                        rd2_wr_file_name,
                        ssake_out_base_name,
                        FRAGMENT_LEN,
                        phred_threshold=20,  # -x
                        n_consec_bases=70,   # -n
                        ascii_offset=33,     # -d
                        min_coverage=5,      # -w
                        n_ovrlap_bases=20,   # -m
                        n_reads_to_call=2,   # -o
                        base_ratio=0.7,      # -r
                        n_bases_to_trim=0,   # -t
                        contig_size=100,     # -z
                        do_track_cvrg=0,     # -c
                        do_ignore_mppng=0,   # -y
                        do_ignore_headr=0,   # -k
                        do_break_ties=0,     # -q
                        do_run_verbose=0)    # -v
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # Assemble paired reads with trusted contigs using SPAdes
        if options.assemble_using_spades:
            start_time = time.time()
            print("Assembling paired reads repeats " +
                  "and with trusted contigs using SPAdes ...", end=" ",
                  flush=True)
            spades_wc_out_dir = BASE_FILE_NAME + "_spades_wc"
            trusted_contigs_fNm = os.path.join(
                ssake_out_base_name, ssake_out_base_name + "_scaffolds.fa")
            if os.path.exists(spades_wc_out_dir):
                shutil.rmtree(spades_wc_out_dir)
            command = utilities.spades(rd1_file_name,
                                       rd2_file_name,
                                       spades_wc_out_dir,
                                       trusted_contigs_fNm=trusted_contigs_fNm)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)

        # Assemble paired reads and with long reads using Unicycler
        if options.assemble_using_unicycler:

            # Read the multi-FASTA SSAKE output file
            seq_rcds = [seq_record for seq_record in SeqIO.parse(
                os.path.join(ssake_out_base_name,
                             ssake_out_base_name + "_scaffolds.fa"), "fasta")]

            # Simulate unpaired reads of the longest sequence`
            start_time = time.time()
            print("Simulating unpaired reads from longest sequence ...",
                  end=" ", flush=True)
            rd_file_name = utilities.simulate_unpaired_reads(
                seq_rcds[0].seq,
                BASE_FILE_NAME,
                number_pairs=NUMBER_READS,
                len_read=LEN_UNPAIRED_READ)
            print("done in {0} s".format(time.time() - start_time), flush=True)

            # Assemble paired reads without repeats and with a long
            # read (SSAKE assembly) using Unicycler
            start_time = time.time()
            print("Assembling paired reads " +
                  "and with a long read using Unicycler ...", end=" ",
                  flush=True)
            unicycler_wc_out_dir = BASE_FILE_NAME + "_unicycler_wc"
            if os.path.exists(unicycler_wc_out_dir):
                shutil.rmtree(unicycler_wc_out_dir)
            command = utilities.unicycler(rd1_file_name,
                                          rd2_file_name,
                                          unicycler_wc_out_dir,
                                          inp_fq_lng_fNm=rd_file_name)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)
