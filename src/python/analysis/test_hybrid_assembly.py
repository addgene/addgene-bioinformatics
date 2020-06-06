from argparse import ArgumentParser
import os
import pickle
import shutil
import time

import utilities


# Set processing parameters
BASE_FILE_NAME = "hybrid"

SEQ_LEN = 10000

K_MER_LEN = 25
EXP_CNT = 16
K_MER_CNT_REP = [EXP_CNT]

NUMBER_PAIRS = 25000
OUTER_DISTANCE = 500
LEN_FIRST_READ = 250
LEN_SECOND_READ = 250

APPEND_DISTANCE = LEN_FIRST_READ

NUMBER_READS = 250
LEN_UNPAIRED_READ = 500

FRAGMENT_LEN = OUTER_DISTANCE


if __name__ == "__main__":

    # Add and parser arguments
    parser = ArgumentParser()
    parser.add_argument('-i', '--initialize',
                        action='store_true',
                        help="create sequence and count k-mers")
    parser.add_argument('-k', '--assemble-using-ssake',
                        action='store_true',
                        help="assemble using SSAKE")
    parser.add_argument('-s', '--assemble-using-spades',
                        action='store_true',
                        help="assemble using SPAdes")
    parser.add_argument('-u', '--assemble-using-unicycler',
                        action='store_true',
                        help="assemble using Unicycler")
    options = parser.parse_args()

    # Create and dump, or load results
    pickle_file_name = BASE_FILE_NAME + ".pickle"
    if not os.path.exists(pickle_file_name) or options.initialize:

        # Create a random sequence with repeats
        start_time = time.time()
        print("Creating random sequence with repeats ...", end=" ", flush=True)
        # Note that repeats need to occur in the interior of the
        # sequence in order to obtain the expected count since we
        # are simapling in a linear sequence. This issue may be
        # resolved with the random rotation added to the paired
        # read simulation
        a_seq = utilities.create_r_seq(
            int((SEQ_LEN - sum(K_MER_CNT_REP) * K_MER_LEN) / 2)
        )
        b_seq, r_k_mers = utilities.create_r_seq_w_rep(
            k_mer_len=K_MER_LEN, k_mer_cnt=K_MER_CNT_REP
        )
        c_seq = utilities.create_r_seq(
            SEQ_LEN - len(a_seq) - len(b_seq)
        )
        d_seq = a_seq[-APPEND_DISTANCE:] + b_seq + c_seq[:APPEND_DISTANCE]
        r_seq = a_seq + b_seq + c_seq
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Random sequence length: {0}".format(len(r_seq)))

        # Simulate paired reads of the entire sequence
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

        # Simulate unpaired reads of the repeat containing sequence
        start_time = time.time()
        print("Simulating unpaired reads from repeat containing sequence ...",
              end=" ", flush=True)
        rd_file_name_q, rd_file_name_a = utilities.simulate_unpaired_reads(
            d_seq,
            BASE_FILE_NAME,
            number_pairs=NUMBER_READS,
            len_read=LEN_UNPAIRED_READ)
        print("done in {0} s".format(time.time() - start_time), flush=True)

        # Create the results dictionary, and dump to a pickle file
        start_time = time.time()
        print("Dumping results ...", end=" ", flush=True)
        results = {}
        results['a_seq'] = a_seq
        results['b_seq'] = b_seq
        results['c_seq'] = c_seq
        results['d_seq'] = d_seq
        results['r_seq'] = r_seq
        results['r_k_mers'] = r_k_mers
        results['rd1_file_name'] = rd1_file_name
        results['rd2_file_name'] = rd2_file_name
        results['rd_file_name_q'] = rd_file_name_q
        results['rd_file_name_a'] = rd_file_name_a
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
        a_seq = results['a_seq']
        b_seq = results['b_seq']
        c_seq = results['c_seq']
        d_seq = results['d_seq']
        r_seq = results['r_seq']
        r_k_mers = results['r_k_mers']
        rd1_file_name = results['rd1_file_name']
        rd2_file_name = results['rd2_file_name']
        rd_file_name_q = results['rd_file_name_q']
        rd_file_name_a = results['rd_file_name_a']
        print("done in {0} s".format(time.time() - start_time), flush=True)

    # Assemble paired reads using SSAKE
    if options.assemble_using_ssake:
        start_time = time.time()
        print("Assembling paired reads using SSAKE ...", end=" ", flush=True)
        ssake_out_base_name = BASE_FILE_NAME + "_ssake"
        if os.path.exists(ssake_out_base_name):
            shutil.rmtree(ssake_out_base_name)
        utilities.ssake(rd1_file_name,
                        rd2_file_name,
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
        print("Assembling using SPAdes ...", end=" ", flush=True)
        spades_wc_out_dir = BASE_FILE_NAME + "_spades"
        if os.path.exists(spades_wc_out_dir):
            shutil.rmtree(spades_wc_out_dir)
        command = utilities.spades(rd1_file_name,
                                   rd2_file_name,
                                   spades_wc_out_dir,
                                   trusted_contigs_fNm=rd_file_name_a)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Command: {0}".format(command), flush=True)

    # Assemble paired reads with long reads using Unicycler
    if options.assemble_using_unicycler:
        start_time = time.time()
        print("Assembling paired reads with long reads using Unicycler ...",
              end=" ", flush=True)
        unicycler_wc_out_dir = BASE_FILE_NAME + "_unicycler_wc"
        if os.path.exists(unicycler_wc_out_dir):
            shutil.rmtree(unicycler_wc_out_dir)
        command = utilities.unicycler(rd1_file_name,
                                      rd2_file_name,
                                      unicycler_wc_out_dir,
                                      inp_fq_lng_fNm=rd_file_name_q)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Command: {0}".format(command), flush=True)
