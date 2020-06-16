from argparse import ArgumentParser
from configparser import ConfigParser
import os
import pickle
import shutil
import time

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

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


"""
Simulates hybrid assembly using simulated paired reads of a known
random sequence, and simulated unpaired reads of a known repetitive,
or random, embedded subsequence.
  - Creates a random sequence, embedding a repetitive or random
    subsequence
  - Simulates paired reads of the random sequence
  - Simulates unpaired reads of the embedded sequence
  - Assembles paired reads with, or without, unpaired reads using
    SPAdes, and aligns
  - Assembles paired reads with, or without, unpaired reads using
    Unicycler, and aligns

"""
if __name__ == "__main__":

    # Add and parser arguments
    parser = ArgumentParser()
    parser.add_argument('case',
                        action='store',
                        help="case identifier suffix")
    parser.add_argument('-i', '--initialize',
                        action='store_true',
                        help="create sequence and simulate reads")
    parser.add_argument('-e', '--embed-repeat',
                        action='store_true',
                        help="embed repeat within linear sequence")
    parser.add_argument('-w', '--use-wgsim',
                        action='store_true',
                        help="use wgsim to simulate reads")
    parser.add_argument('-s', '--assemble-using-spades',
                        action='store_true',
                        help="assemble using SPAdes")
    parser.add_argument('-u', '--assemble-using-unicycler',
                        action='store_true',
                        help="assemble using Unicycler")
    parser.add_argument('-l', '--include-long-reads',
                        action='store_true',
                        help="include long reads during asssembly")
    parser.add_argument('-f', '--force-assembly',
                        action='store_true',
                        help="run assemblers when output directory exists")
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
    case_file_name = BASE_FILE_NAME + "_" + options.case
    pickle_file_name = case_file_name + ".pickle"
    if not os.path.exists(pickle_file_name) or options.initialize:

        # Create a random sequence
        start_time = time.time()
        print("Creating random sequence ...", end=" ", flush=True)
        rep_len = sum(K_MER_CNT_REP) * K_MER_LEN
        a_seq = utilities.create_r_seq(
            int((SEQ_LEN - rep_len) / 2)
        )
        if options.embed_repeat:
            b_seq, r_k_mers = utilities.create_r_seq_w_rep(
                k_mer_len=K_MER_LEN, k_mer_cnt=K_MER_CNT_REP
            )
        else:
            b_seq = utilities.create_r_seq(
                rep_len
            )
            r_k_mers = []
        c_seq = utilities.create_r_seq(
            SEQ_LEN - len(a_seq) - len(b_seq)
        )
        d_seq = a_seq[-APPEND_DISTANCE:] + b_seq + c_seq[:APPEND_DISTANCE]
        r_seq = a_seq + b_seq + c_seq
        r_seq_rcd = SeqRecord(
            r_seq, id="1", name="reference", description="complete plasmid")
        r_seq_fNm = case_file_name + ".fasta"
        SeqIO.write(r_seq_rcd, r_seq_fNm, "fasta")
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Random sequence length: {0}".format(len(r_seq)))

        # Simulate paired reads of the random sequence
        start_time = time.time()
        if options.use_wgsim:
            print(
                "Simulating paired reads from random sequence using wgsim ...",
                end=" ", flush=True)
            rd1_file_name, rd2_file_name = utilities.simulate_paired_reads_wgsim(
                r_seq,
                case_file_name,
                number_pairs=NUMBER_PAIRS,
                len_first_read=LEN_FIRST_READ,
                len_second_read=LEN_SECOND_READ,
                outer_distance=OUTER_DISTANCE)
        else:
            print(
                "Simulating paired reads from random sequence ...",
                end=" ", flush=True)
            rd1_file_name, rd2_file_name = utilities.simulate_paired_reads(
                r_seq,
                case_file_name,
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
            case_file_name,
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

    # Assemble paired reads using SPAdes
    spades_out_dir = case_file_name + "_spades"
    if (options.assemble_using_spades and (
        not os.path.exists(spades_out_dir) or options.force_assembly
    )):
        start_time = time.time()
        print("Assembling using SPAdes ...", end=" ", flush=True)
        if os.path.exists(spades_out_dir):
            shutil.rmtree(spades_out_dir)
        if options.include_long_reads:
            command = utilities.spades(rd1_file_name,
                                       rd2_file_name,
                                       spades_out_dir,
                                       trusted_contigs_fNm=rd_file_name_a)
        else:
            command = utilities.spades(rd1_file_name,
                                       rd2_file_name,
                                       spades_out_dir)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Command: {0}".format(command), flush=True)

    # Read the multi-FASTA SPAdes output file
    if os.path.exists(spades_out_dir):
        spades_seq_rcds = [seq_record for seq_record in SeqIO.parse(
            os.path.join(spades_out_dir,
                         "contigs.fasta"), "fasta")]
        if len(spades_seq_rcds) > 0:
            print("SPAdes score: {0}".format(
                aligner.score(r_seq + r_seq, spades_seq_rcds[0].seq)))

    # Assemble paired reads using Unicycler
    unicycler_out_dir = case_file_name + "_unicycler"
    if (options.assemble_using_unicycler and (
        not os.path.exists(unicycler_out_dir) or options.force_assembly
    )):
        start_time = time.time()
        print("Assembling using Unicycler ...", end=" ", flush=True)
        if os.path.exists(unicycler_out_dir):
            shutil.rmtree(unicycler_out_dir)
        if options.include_long_reads:
            command = utilities.unicycler(rd1_file_name,
                                          rd2_file_name,
                                          unicycler_out_dir,
                                          inp_fq_lng_fNm=rd_file_name_q)
        else:
            command = utilities.unicycler(rd1_file_name,
                                          rd2_file_name,
                                          unicycler_out_dir)
        print("done in {0} s".format(time.time() - start_time), flush=True)
        print("Command: {0}".format(command), flush=True)

    # Read the multi-FASTA Unicycler output file
    if os.path.exists(unicycler_out_dir):
        unicycler_seq_rcds = [seq_record for seq_record in SeqIO.parse(
            os.path.join(unicycler_out_dir,
                         "assembly.fasta"), "fasta")]
        if len(unicycler_seq_rcds) > 0:
            print("Unicycler score: {0}".format(
                aligner.score(r_seq + r_seq, unicycler_seq_rcds[0].seq)))
