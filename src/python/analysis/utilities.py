import gzip
import logging
import math
import os
import random
import re

from Bio import Align
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import docker
import numpy as np

logger = logging.getLogger(__name__)


def create_r_seq_str(seq_len):
    """Create a string of randomly selected nucleotides.

    Parameters
    ----------
    seq_len : int
        the length of the string to creeate

    Returns
    -------
    str
        the sequence string
    """
    r_seq_str = ""
    for count in range(seq_len):
        r_seq_str += random.choice("CGTA")
    return r_seq_str


def create_r_seq(seq_len):
    """Creates a sequence of randomly selected nucleotides.

    Parameters
    ----------
    seq_len : int
        the length of the sequence to create

    Returns
    -------
    Bio.Seq.Seq
        the sequence
    """
    r_seq_str = create_r_seq_str(seq_len)
    r_seq = Seq(r_seq_str, IUPAC.unambiguous_dna)
    return r_seq


def create_r_seq_w_rep(k_mer_len=25, k_mer_cnt=[2 ** (n + 1) for n in range(8)]):
    """Creates a k-mer set of randomly selected nucleotides and uses
    the set to construct a sequence with the specified number of
    repeats.

    Parameters
    ----------
    k_mer_len : int
        the length of the k_mer to create
    k_mer_cnt : [int]
        the number of times to repeat each k_mer

    Returns
    -------
    (Bio.Seq.Seq, [Bio.Seq.Seq])
        the sequence and list of repeated k-mers
    """
    r_k_mers = []
    r_seq = create_r_seq(0)
    for kMC in k_mer_cnt:
        r_k_mer = create_r_seq(k_mer_len)
        r_k_mers.append(r_k_mer)
        for iKM in range(kMC):
            r_seq += r_k_mer
    return r_seq, r_k_mers


# TODO: Add error rate, outer distance standard deviation, indel
# fraction, indel extension probablility, and random seed:
# error_rate=0.020, standard_deviation=50, indel_fraction=0.15,
# indel_extended_prob=0.30, and random_seed=0, then move to utilities
def simulate_paired_reads(
    seq,
    seq_nm,
    number_pairs=25000,
    len_first_read=250,
    len_second_read=250,
    outer_distance=500,
):
    rd1_seq_rcds = []
    rd2_seq_rcds = []
    seq_len = len(seq)
    for iP in range(number_pairs):

        # Create a random rotation of the sequence
        i_seq = random.randint(0, seq_len)
        r_seq = seq[seq_len - i_seq : seq_len] + seq[0 : seq_len - i_seq]

        # Create the first read starting from a random index, and
        # setting random PHRED quality scores
        i_rd1 = random.randint(0, seq_len - outer_distance - 1)
        rd1_seq_rcd = SeqRecord(
            r_seq[i_rd1 : i_rd1 + len_first_read],
            id=str(iP),
            name="one",
            description="first read of read pair",
        )
        rd1_seq_rcd.letter_annotations["phred_quality"] = (
            np.random.randint(20, 50, len_first_read)
        ).tolist()
        rd1_seq_rcds.append(rd1_seq_rcd)

        # Create the second read starting from an index giving the
        # correct outer distance, and setting random PHRED quality
        # scores
        i_rd2 = i_rd1 + outer_distance - len_second_read
        rd2_seq_rcd = SeqRecord(
            r_seq[i_rd2 + len_second_read : i_rd2 : -1],
            id=str(iP),
            name="two",
            description="second read of read pair",
        )
        rd2_seq_rcd.letter_annotations["phred_quality"] = (
            np.random.randint(20, 50, len_second_read)
        ).tolist()
        rd2_seq_rcds.append(rd2_seq_rcd)

    rd1_fNm = seq_nm + "_rd1.fastq"
    rd2_fNm = seq_nm + "_rd2.fastq"
    SeqIO.write(rd1_seq_rcds, rd1_fNm, "fastq")
    SeqIO.write(rd2_seq_rcds, rd2_fNm, "fastq")
    return rd1_fNm, rd2_fNm


def simulate_unpaired_reads(seq, seq_nm, number_pairs=25000, len_read=500):
    rd_seq_rcds = []
    seq_len = len(seq)
    for iP in range(number_pairs):

        # Create a random rotation of the sequence
        i_seq = random.randint(0, seq_len)
        r_seq = seq[seq_len - i_seq : seq_len] + seq[0 : seq_len - i_seq]

        # Create the read starting from a random index, and
        # setting random PHRED quality scores
        i_rd = random.randint(0, seq_len - len_read - 1)
        rd_seq_rcd = SeqRecord(
            r_seq[i_rd : i_rd + len_read],
            id=str(iP),
            name="unpaired",
            description="unpaired read",
        )
        rd_seq_rcd.letter_annotations["phred_quality"] = (
            np.random.randint(20, 50, len_read)
        ).tolist()
        rd_seq_rcds.append(rd_seq_rcd)

    rd_fNm = seq_nm + "_rd.fastq"
    SeqIO.write(rd_seq_rcds, rd_fNm, "fastq")
    return rd_fNm


def create_aligner(config):
    """Creates a pairwise aligner with the specified configuration.

    Parameters
    ----------
    config : dct
        pairwise aligner configuration
    Returns
    -------
    Align.PairwiseAligner
        the Biopython pairwise aligner
    """
    # Create a pairwise aligner with the specified configuration
    aligner = Align.PairwiseAligner()
    aligner.match_score = float(config["match_score"])
    aligner.mismatch_score = float(config["mismatch_score"])
    aligner.target_open_gap_score = float(config["target_open_gap_score"])
    aligner.target_extend_gap_score = float(config["target_extend_gap_score"])
    aligner.target_left_open_gap_score = float(config["target_left_open_gap_score"])
    aligner.target_left_extend_gap_score = float(config["target_left_extend_gap_score"])
    aligner.target_right_open_gap_score = float(config["target_right_open_gap_score"])
    aligner.target_right_extend_gap_score = float(
        config["target_right_extend_gap_score"]
    )
    aligner.query_open_gap_score = float(config["query_open_gap_score"])
    aligner.query_extend_gap_score = float(config["query_extend_gap_score"])
    aligner.query_left_open_gap_score = float(config["query_left_open_gap_score"])
    aligner.query_left_extend_gap_score = float(config["query_left_extend_gap_score"])
    aligner.query_right_open_gap_score = float(config["query_right_open_gap_score"])
    aligner.query_right_extend_gap_score = float(config["query_right_extend_gap_score"])
    aligner.mode = config["mode"]

    return aligner


def count_k_mers_in_seq(seq, seq_id=0, k_mer_len=25, k_mers=None):
    """Count k-mers of a specified length in a sequence with a
    specified identifier.

    Parameters
    ----------
    seq : Bio.Seq.Seq
        the sequence in which to count
    seq_id : int
        the sequence identifer
    k_mer_len : int
        the length of k-mers to count
    k_mers : None or dct
        None or dictionary returned by this method

    Returns
    -------
    dct
        dictionay containing k-mer keys and counts and source sequence
        identifiers values
    """
    if k_mers is None:
        k_mers = {}
    seq_len = len(seq)
    for i_seq in range(seq_len - k_mer_len + 1):
        k_mer = str(seq[i_seq : i_seq + k_mer_len])
        if k_mer not in k_mers:
            k_mers[k_mer] = {}
            k_mers[k_mer]["src"] = set([seq_id])
            k_mers[k_mer]["cnt"] = 1
        else:
            k_mers[k_mer]["src"].add(seq_id)
            k_mers[k_mer]["cnt"] += 1
    return k_mers


def count_k_mers_in_rds(rd_fNm, k_mer_len=25, k_mers=None, seq_rcds=None):
    """Count k-mers of a specified length in a gzipped file of reads
    in the specified format.

    Parameters
    ----------
    rd_fNm : str
        the name of the gzipped read file
    k_mer_len : int
        the length of k-mers to count
    k_mers : None or dct
        None or dictionary returned by count_k_mers_in_seq()
    seq_rcds : None or lst
        None or list of sequence records in which k-mers were counted

    Returns
    -------
    (dct, [Bio.SeqRecord.SeqRecord])
        dictionay containing k-mer keys and count and source sequence
        record index values, and list of sequence records in which
        k-mers were counted
    """
    if k_mers is None:
        k_mers = {}
    if seq_rcds is None:
        seq_rcds = []
    read_format, is_gzipped = get_bio_read_format(rd_fNm)
    _open = open
    if is_gzipped:
        _open = gzip.open
    with _open(rd_fNm, "rt") as f:
        seq_rcd_gen = SeqIO.parse(f, format=read_format)
        i_seq = -1
        for seq_rcd in seq_rcd_gen:
            i_seq += 1
            k_mers = count_k_mers_in_seq(seq_rcd.seq, i_seq, k_mers=k_mers)
            seq_rcds.append(seq_rcd)
    return k_mers, seq_rcds


def write_k_mer_counts_in_rds(k_mers_in_rd1, k_mers_in_rd2, k_mer_counts_fNm):
    """Write the k-mer counts corresponding to paired reads.

    Parameters
    ----------
    k_mers_in_rd1 : dct
        dictionary returned by count_k_mers_in_seq()
    k_mers_in_rd2 : dct
        dictionary returned by count_k_mers_in_seq()
    k_mer_counts_fNm
        file name to which to write k-mers and their counts
    """
    with open(k_mer_counts_fNm, "w") as f:
        k_mer_set_rd1 = set(k_mers_in_rd1.keys())
        k_mer_set_rd2 = set(k_mers_in_rd2.keys())
        k_mer_set = k_mer_set_rd1.union(k_mer_set_rd2)
        for k_mer in sorted(k_mer_set):
            k_mer_cnt_rd1 = 0
            if k_mer in k_mer_set_rd1:
                k_mer_cnt_rd1 = k_mers_in_rd1[k_mer]["cnt"]
            k_mer_cnt_rd2 = 0
            if k_mer in k_mer_set_rd2:
                k_mer_cnt_rd2 = k_mers_in_rd2[k_mer]["cnt"]
            f.write(
                "{0} {1:10d} {1:10d}\n".format(
                    k_mer, int(k_mer_cnt_rd1), int(k_mer_cnt_rd2)
                )
            )


def read_k_mer_counts(k_mer_counts_fNm, seq_id=0, k_mers=None):
    """Read k-mer counts from the specified file.

    Parameters
    ----------
    k_mers_counts_fNm : str
        the name of the file containing k-mers and their counts
    seq_id : int
        the sequence identifer
    k_mers : None or dct
        None or dictionary returned by count_k_mers_in_seq()

    Returns
    -------
    dct
        dictionay containing k-mer keys, and counts and source sequence
        identifiers values
    """
    if k_mers is None:
        k_mers = {}
    with open(k_mer_counts_fNm) as f:
        for ln in f:
            flds = ln.split()
            k_mer = flds[0]
            cnt_rd1 = int(flds[1])
            cnt_rd2 = int(flds[2])
            if k_mer in k_mers:
                print("Second occurance of k-mer: {0} unexpected".format(k_mer))
            else:
                k_mers[k_mer] = {}
                k_mers[k_mer]["src"] = set([seq_id])
                k_mers[k_mer]["rd1"] = cnt_rd1
                k_mers[k_mer]["rd2"] = cnt_rd2
    return k_mers


def rotate_seqs(a_seq, o_seq):
    """Finds the cyclic rotation of a_seq (or an approximation of it)
    that minimizes the blockwise q-gram distance from o_seq using
    Docker image "ralatsdio/csc:v0.1.0".

    Parameters
    ----------
    a_seq : Bio.Seq.Seq
        thought of as the reference sequence
    o_seq : Bio.Seq.Seq
        thought of as the candidate sequence

    Returns
    -------
    tpl
        tuple containing the rotated reference and candidate sequences
    """
    # Define image, Docker run parameters, and CSC command
    image = "ralatsdio/csc:v0.1.0"
    # [-a] `DNA' or `RNA' for nucleotide sequences or `PROT' for
    # protein sequences
    alphabet = "DNA"
    # [-m] `hCSC' for heuristic, `nCSC' for naive and `saCSC' for
    # suffix-array algorithm
    method = "saCSC"
    # (Multi)FASTA input filename.
    input_file = "csc_input.fasta"
    # Output filename for the rotated sequences.
    output_file = "csc_output.fasta"
    # [-l] The length of each block
    block_length = str(math.ceil(math.sqrt(len(a_seq))))
    # [-q] The q-gram length
    q_length = str(math.ceil(math.log(len(a_seq)) / math.log(4)))
    # [-P] The number of blocks of length l to use to refine the
    # results of saCSC by (e.g. 1.0)
    blocks_refine = "1.0"
    # [-O] The gap open penalty is the score taken away when a gap is
    # created. The best value depends on the choice of comparison
    # matrix.  The default value assumes you are using the EBLOSUM62
    # matrix for protein sequences, and the EDNAFULL matrix for
    # nucleotide sequences. Floating point number from 1.0 to
    # 100.0. (default: 10.0)
    gap_open_penalty = "10.0"
    # [-E] The gap extension penalty is added to the standard gap
    # penalty for each base or residue in the gap. This is how long
    # gaps are penalized. Floating point number from 0.0 to 10.0.
    # (default: 0.5)
    gap_extend_penalty = "0.5"
    command = " ".join(
        [
            "csc",
            "-m",
            method,
            "-a",
            alphabet,
            "-i",
            input_file,
            "-o",
            output_file,
            "-q",
            q_length,
            "-l",
            block_length,
            "-P",
            blocks_refine,
            "-O",
            gap_open_penalty,
            "-E",
            gap_extend_penalty,
        ]
    )
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Write the multi-FASTA input file
    SeqIO.write(
        [
            SeqRecord(a_seq, id="id_a", name="name_a", description="reference"),
            SeqRecord(o_seq, id="id_o", name="name_o", description="offset"),
        ],
        os.path.join(hosting_dir, input_file),
        "fasta",
    )

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )

    # Read the multi-FASTA output file
    seq_rcds = [
        seq_record
        for seq_record in SeqIO.parse(os.path.join(hosting_dir, output_file), "fasta")
    ]

    # Assign and return the rotated sequences
    r_a_seq = seq_rcds[0].seq
    r_o_seq = seq_rcds[1].seq

    return (r_a_seq, r_o_seq)


# TODO: Complete in order to validate use of wgsim
def simulate_paired_reads_wgsim(
    seq,
    seq_nm,
    number_pairs=25000,
    len_first_read=250,
    len_second_read=250,
    outer_distance=500,
):
    seq_rcd = SeqRecord(seq, id="0", name="base", description="reference")
    seq_fNm = seq_nm + ".fasta"
    rd1_fNm = seq_nm + "_rd1.fastq"
    rd2_fNm = seq_nm + "_rd2.fastq"
    SeqIO.write(seq_rcd, seq_fNm, "fasta")
    wgsim(
        seq_fNm,
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
        haplotype_mode=True,
    )
    return rd1_fNm, rd2_fNm


def kmc(
    read_file_names,
    database_file_name,
    k_mer_length=25,
    signature_length=9,
    count_min=2,
    max_count=255,
    count_max=1e9,
    canonical_form=True,
):
    """
    Counts k-mers.

    Documentation for kmc ver. 3.1.1 (2019-05-19):

    Usage:
        kmc [options] <input_file_name> <output_file_name> <working_directory>
        kmc [options] <@input_file_names> <output_file_name> <working_directory>

    Options:
        -v - verbose mode (shows all parameter settings); default: false
        -k<len> - k-mer length (k from 1 to 256; default: 25)
        -m<size> - max amount of RAM in GB (from 1 to 1024); default: 12
        -sm - use strict memory mode (memory limit from -m<n> switch
            will not be exceeded)
        -p<par> - signature length (5, 6, 7, 8, 9, 10, 11); default: 9
        -f<a/q/m/bam> - input in FASTA format (-fa), FASTQ format
            (-fq), multi FASTA (-fm) or BAM (-fbam); default: FASTQ
        -ci<value> - exclude k-mers occurring less than <value> times
            (default: 2)
        -cs<value> - maximal value of a counter (default: 255)
        -cx<value> - exclude k-mers occurring more of than <value>
            times (default: 1e9)
        -b - turn off transformation of k-mers into canonical form
        -r - turn on RAM-only mode
        -n<value> - number of bins
        -t<value> - total number of threads (default: no. of CPU cores)
        -sf<value> - number of FASTQ reading threads
        -sp<value> - number of splitting threads
        -sr<value> - number of threads for 2nd stage
        -j<file_name> - file name with execution summary in JSON format
        -w - without output

    Parameters:
        input_file_name - single file in specified (-f switch) format
            (gziped or not)
        @input_file_names - file name with list of input files in
            specified (-f switch) format (gziped or not)
        output_file_name - path to output file
        working_directory - directory containing intermediate files

    Examples:
        $ kmc -k27 -m24 NA19238.fastq NA.res /data/kmc_tmp_dir/
        $ kmc -k27 -m24 @files.lst NA.res /data/kmc_tmp_dir/
    """
    # Determine the read file format
    read_format = get_kmc_read_format(read_file_names)

    # Determine whether to submit read file names as a string, or file
    if len(read_file_names) > 1:
        with open("kmc_input.txt", "w") as f:
            for input_file_name in read_file_names:
                f.write("{0}\n".format(input_file_name))
        input_str = "@kmc_input.txt"
    else:
        input_str = read_file_names[0]

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define kmc command
    command = " ".join(
        [
            "kmc",
            "-k" + str(k_mer_length),
            "-p" + str(signature_length),
            "-ci" + str(int(count_min)),
            "-cs" + str(int(max_count)),
            "-cx" + str(int(count_max)),
        ]
    )
    if not canonical_form:
        command = " ".join([command, "-b"])
    command = " ".join(
        [command, read_format, input_str, database_file_name, working_dir,]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )


def kmc_transform(
    inp_database_file_name,
    operation,
    out_database_file_name,
    inp_count_min=2,
    inp_count_max=1e9,
    out_count_min=2,
    out_count_max=1e9,
    is_sorted=False,
):
    """
    Transforms single input database to output (text file or KMC database).

    Documentation for kmc_tools transform ver. 3.1.1 (2019-05-19):

    Usage:
        kmc_tools [options] transform <input> [input_params] \
            <oper1 [oper_params1] output1 [output_params1]> \
            [<oper2 [oper_params2] output2 [output_params2]> ...]

    Options:
        -t<value> - total number of threads (default: no. of CPU cores)
        -v - enable verbose mode (shows some information) (default: false)
        -hp - hide percentage progress (default: false)

    Transform parameters:
        input - path to database generated by KMC
        oper1, oper2, ..., operN - transform operation name
        output1, output2, ..., outputN - paths to output

    Input parameters:
        -ci<value> - exclude k-mers occurring less than <value> times
        -cx<value> - exclude k-mers occurring more of than <value> times

    Transform operations:
        sort - converts database produced by KMC2.x to KMC1.x database
            format (which contains k-mers in sorted order)
        reduce - exclude too rare and too frequent k-mers
        compact - remove counters of k-mers
        histogram - produce histogram of k-mers occurrences
        dump - produce text dump of kmc database
        set_counts <value> - set all k-mer counts to specific value

    Sort and Reduce output parameters:
        -ci<value> - exclude k-mers occurring less than <value> times
        -cx<value> - exclude k-mers occurring more of than <value> times
        -cs<value> - maximal value of a counter

    Histogram output parameters:
        -ci<value> - minimum value of counter to be stored in the otput file
        -cx<value> - maximum value of counter to be stored in the otput file

    Dump operation parameters:
        -s - sorted output

    Example:
        kmc_tools transform db \
            reduce err_kmers -cx10 \
            reduce valid_kmers -ci11 \
            histogram histo.txt \
            dump dump.txt
    """
    # Check input arguments
    if operation not in ["reduce", "histogram", "dump"]:
        raise (Exception("Operation is not implemented"))

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define kmc_tools command
    command = " ".join(
        [
            "kmc_tools",
            "transform",
            inp_database_file_name,
            "-ci" + str(int(inp_count_min)),
            "-cx" + str(int(inp_count_max)),
            operation,
            out_database_file_name,
            "-ci" + str(int(out_count_min)),
            "-cx" + str(int(out_count_max)),
        ]
    )
    if operation == "dump" and is_sorted:
        command = " ".join([command, "-s"])

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )


# TODO: Validate
def kmc_simple(
    inp_database_file_name_a,
    operation,
    inp_database_file_name_b,
    out_database_file_name,
    inp_count_min_a=2,
    inp_count_max_a=1e9,
    inp_count_min_b=2,
    inp_count_max_b=1e9,
    out_count_min=2,
    out_max_count=255,
    out_count_max=1e9,
    out_calc_mode="",
):
    """
    Performs set operation on two input KMC databases.

    Documentation for kmc_tools simple ver. 3.1.1 (2019-05-19):

    Usage:
        kmc_tools [options] simple <input1 [input1_params]> <input2 [input2_params]> \
            <oper1 output1 [output_params1]> \
            [<oper2 output2 [output_params2]> ...]

    Options:
        -t<value> - total number of threads (default: no. of CPU cores)
        -v - enable verbose mode (shows some information) (default: false)
        -hp - hide percentage progress (default: false)

    Simple parameters:
        input1, input2 - paths to databases generated by KMC
        oper1, oper2, ..., operN - set operations to be performed on
            input1 and input2
        output1, output2, ..., outputN - paths to output k-mer databases

    Input parameters:
        -ci<value> - exclude k-mers occurring less than <value> times
        -cx<value> - exclude k-mers occurring more of than <value> times

    Simple operations:
        intersect - output database will contains only k-mers that are
            present in both input sets
        union - output database will contains each k-mer present in
            any of input sets
        kmers_subtract - difference of input sets based on
            k-mers. Output database will contains only k-mers that are
            present in first input set but absent in the second one
        counters_subtract - difference of input sets based on k-mers
            and their counters (weaker version of
            kmers_subtract). Output database will contains all k-mers
            that are present in first input, beyond those for which
            counter operation will lead to remove (i.e. counter equal
            to 0 or negative number)
        reverse_kmers_subtract - same as kmers_subtract but treat
            input2 as first and input1 as second
        reverse_counters_subtract - same as counters_subtract but
            treat input2 as first and input1 as second

    Output parameters:
        -ci<value> - exclude k-mers occurring less than <value> times
        -cx<value> - exclude k-mers occurring more of than <value> times
        -cs<value> - maximal value of a counter
        -oc<value> - redefine counter calculation mode for equal k-mers

    Available values:
        min - get lower value of a k-mer counter (default value for
            intersect operation)
        max - get upper value of a k-mer counter
        sum - get sum of counters from both databases
        diff - get difference between counters (default for
            counters_subtract operation)
        left - get counter from first database (input1)
        right - get counter from second database (input2)

    Example:
        kmc_tools simple kmers1 -ci3 -cx70000 kmers2 \
            union kmers1_kmers2_union -cs65536 -ocfirst \
            intersect intersect_kmers1_kmers2 \
            intersect intersect_max_kmers1_kmers2 -ocmax

    kmc_tools ver. 3.1.1 (2019-05-19)
    """
    # Check input arguments
    if operation not in [
        "intersect",
        "union",
        "kmers_subtract",
        "counters_subtract",
        "reverse_kmers_subtract",
        "reverse_counters_subtract",
    ]:
        raise (Exception("Invalid operation: {0}".format(operation)))
    if out_calc_mode not in ["", "min", "max", "sum", "diff", "left", "right"]:
        raise (Exception("Invalid output calculation mode: {0}".format(out_calc_mode)))

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define kmc_tools command
    command = " ".join(
        [
            "kmc_tools",
            "simple",
            inp_database_file_name_a,
            "-ci" + str(int(inp_count_min_a)),
            "-cx" + str(int(inp_count_max_a)),
            inp_database_file_name_b,
            "-ci" + str(int(inp_count_min_b)),
            "-cx" + str(int(inp_count_max_b)),
            operation,
            out_database_file_name,
            "-ci" + str(int(out_count_min)),
            "-cs" + str(int(out_max_count)),
            "-cx" + str(int(out_count_max)),
        ]
    )
    if out_calc_mode != "":
        command = " ".join([command, "-oc{0}".format(out_calc_mode)])

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )


# TODO: Validate
def kmc_complex():
    """
    Performs set operations for more than two input k-mers sets.

    Documentation for kmc_tools complex ver. 3.1.1 (2019-05-19):

    Usage:
        kmc_tools [options] complex <operations_definition_file>

    Options:
        -t<value> - total number of threads (default: no. of CPU cores)
        -v - enable verbose mode (shows some information) (default: false)
        -hp - hide percentage progress (default: false)

    Complex parameters:
        operations_definition_file - path to file which define input
            sets and operations.

    Operations definition file syntax:

        INPUT:
        <input1>=<input1_db_path> [params]
        <input2>=<input2_db_path> [params]
        ...
        <inputN>=<inputN_db_path> [params]
        OUTPUT:
        <out_db_path>=<ref_input><oper[c_mode]><ref_input>[<oper[c_mode]><ref_input>[...]
        [OUTPUT_PARAMS:
        <output_params>]

        Input:
            input1, input2, ..., inputN - names of inputs used to
                define equation
            input1_db, input2_db_path, ..., inputN_db_path - paths to
                k-mers sets

        Input parameters:
            -ci<value> - exclude k-mers occurring less than <value> times
            -cx<value> - exclude k-mers occurring more of than <value> times

        Output:
            out_db_path - path to output database
            ref_input - one of input1, input2, ..., inputN
            oper - one of {*,-,~,+}, which refers to {intersect,
                kmers_subtract, counters_subtract, union}. Operator *
                has the highest priority. Other operators has equal
                priorities. Order of operations can be changed with
                parentheses. For {*,~,+} it is possible to redefine
                counter calculation mode ([c_mode]). Available values:
                min, max, diff, sum, left, right (detailet description
                available in simple help message).

        Output parameters:
            -ci<value> - exclude k-mers occurring less than <value> times
            -cx<value> - exclude k-mers occurring more of than <value> times
            -cs<value> - maximal value of a counter

        Example:
            INPUT:
            set1 = kmc_o1 -ci5
            set2 = kmc_o2
            set3 = kmc_o3 -ci10 -cx100
            OUTPUT:
            result = (set3 + min set1) * right set2

    kmc_tools ver. 3.1.1 (2019-05-19)
    """
    raise (NotImplementedError("utility.kmc_complex() has not been implemented"))


# TODO: Validate
def kmc_filter(
    inp_database_file_name,
    inp_read_file_name,
    out_read_file_name,
    trim_reads=False,
    hard_mask=False,
    inp_db_count_min=2,
    inp_db_count_max=1e9,
    inp_rd_count_min=2,
    inp_rd_count_max=1e9,
):
    """
    Filters out reads with too small number of k-mers.

    Documentation for kmc_tools filter ver. 3.1.1 (2019-05-19):

    Usage:
        kmc_tools [options] filter [filter_options] <kmc_input_db> [kmc_input_db_params] \
            <input_read_set> [input_read_set_params] \
            <output_read_set> [output_read_set_params]
    Options:
        -t<value> - total number of threads (default: no. of CPU cores)
        -v - enable verbose mode (shows some information) (default: false)
        -hp - hide percentage progress (default: false)

    Filter options:
        -t - trim reads on first invalid k-mer instead of remove entirely
        -hm - hard mask invalid k-mers in a read

    Filter parameters:
        kmc_input_db - path to database generated by KMC
        input_read_set - path to input set of reads
        output_read_set - path to output set of reads

    Database parameters:
        -ci<value> - exclude k-mers occurring less than <value> times
        -cx<value> - exclude k-mers occurring more than <value> times

    Input read set parameters:
        -ci<value> - remove reads containing less k-mers than
            value. It can be integer or floating number in range
            [0.0;1.0] (default: 2)
        -cx<value> - remove reads containing more k-mers than
            value. It can be integer or floating number in range
            [0.0;1.0] (default: 1e9)
        -f<a/q> - input in FASTA format (-fa), FASTQ format (-fq);
            default: FASTQ

    Output read set parameters:
        -f<a/q> - output in FASTA format (-fa), FASTQ format (-fq);
            default: same as input

    Example:
        kmc_tools filter kmc_db -ci3 input.fastq -ci0.5 -cx1.0 filtered.fastq
        kmc_tools filter kmc_db input.fastq -ci10 -cx100 filtered.fastq
    """
    # Determine the input file format
    read_format = get_kmc_read_format([inp_read_file_name])

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define kmc command
    command = " ".join(["kmc_tools", "filter"])
    # TODO: Determine if these options are mutually exclusive
    if trim_reads:
        command = " ".join([command, "-t"])
    if hard_mask:
        command = " ".join([command, "-hm"])
    command = " ".join(
        [
            command,
            inp_database_file_name,
            "-ci" + str(int(inp_db_count_min)),
            "-cx" + str(int(inp_db_count_max)),
            inp_read_file_name,
            "-ci" + str(int(inp_rd_count_min)),
            "-cx" + str(int(inp_rd_count_max)),
            read_format,
            out_read_file_name,
            read_format,
        ]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )


def get_kmc_read_format(read_file_names):
    """
    Determine a valid read format implied by the read file name
    extension(s).

    Parameters
    ----------
    read_file_names : [str]
        list of read file names

    Returns
    -------
    str
        valid read format: FASTA is "-fa", FASTQ is "-fq", multi
        FASTA is "-fm", or BAM is "-fbam"; default is "-fq"
    """
    valid_read_format = ""
    p_fa = re.compile(r"\.f(ast)?a(\.gz)?$")
    p_fq = re.compile(r"\.f(ast)?q(\.gz)?$")
    p_bm = re.compile(r"\.bam$")
    p_ss = re.compile(r"^>")
    for input_file_name in read_file_names:
        if p_fa.search(input_file_name) is not None:  # FASTA

            # Count number of sequences
            with open(input_file_name, "r") as f:
                n_seq = 0
                for ln in f:
                    if p_ss.search(ln) is not None:
                        n_seq += 1
            if n_seq == 0:
                raise (Exception("No sequence found in FASTA file"))

            elif n_seq == 1:  # FASTA
                read_format = "-fa"

            else:  # n_seq > 1 -- multi FASTA
                read_format = "-fm"

        elif p_fq.search(input_file_name) is not None:
            read_format = "-fq"  # FASTQ

        elif p_bm.search(input_file_name) is not None:
            read_format = "-fbam"  # BAM

        else:
            raise (Exception("Unknown sequence file format"))

        if valid_read_format == "":
            valid_read_format = read_format

        elif read_format != valid_read_format:
            raise (Exception("Input files must use the same format"))

    return valid_read_format


def get_bio_read_format(read_file_name):
    """
    Determine a valid read format implied by the read file name
    extension.

    Parameters
    ----------
    read_file_name : str
        list of read file names

    Returns
    -------
    (str, boolean)
        valid read format: "fasta" or "fastq", and flag indicating
        gzip extension, or not
    """
    read_format = ""
    is_gzipped = False
    p_fa = re.compile(r"\.f(ast)?a(\.gz)?$")
    p_fq = re.compile(r"\.f(ast)?q(\.gz)?$")
    p_gz = re.compile(r"(\.gz)$")
    if p_fa.search(read_file_name) is not None:
        read_format = "fasta"

    elif p_fq.search(read_file_name) is not None:
        read_format = "fastq"

    else:
        raise (Exception("Unknown read file format"))
    if p_gz.search(read_file_name) is not None:
        is_gzipped = True
    return read_format, is_gzipped


def wgsim(
    inp_fa_fNm,
    out_fq_one_fNm,
    out_fq_two_fNm,
    error_rate=0.020,
    outer_distance=500,
    standard_deviation=50,
    number_pairs=1000000,
    len_first_read=70,
    len_second_read=70,
    mutation_rate=0.0010,
    indel_fraction=0.15,
    indel_extended_prob=0.30,
    random_seed=0,
    ambiguous_base_frac=0.05,
    haplotype_mode=False,
):
    """
    Simulating sequence reads from a reference genome.

    Documentation for samtools ver. 1.10

    Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>

    Options: -e FLOAT      base error rate [0.020]
             -d INT        outer distance between the two ends [500]
             -s INT        standard deviation [50]
             -N INT        number of read pairs [1000000]
             -1 INT        length of the first read [70]
             -2 INT        length of the second read [70]
             -r FLOAT      rate of mutations [0.0010]
             -R FLOAT      fraction of indels [0.15]
             -X FLOAT      probability an indel is extended [0.30]
             -S INT        seed for random generator [0, use the current time]
             -A FLOAT      discard if the fraction of ambiguous bases higher than FLOAT [0.05]
             -h            haplotype mode
    """
    # Define image, and Docker run parameters
    image = "ralatsdio/samtools:v1.10"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define kmc command
    command = " ".join(
        [
            "wgsim",
            "-e" + str(error_rate),
            "-d" + str(outer_distance),
            "-s" + str(standard_deviation),
            "-N" + str(number_pairs),
            "-1" + str(len_first_read),
            "-2" + str(len_second_read),
            "-r" + str(mutation_rate),
            "-R" + str(indel_fraction),
            "-X" + str(indel_extended_prob),
            "-S" + str(random_seed),
            "-A" + str(ambiguous_base_frac),
        ]
    )
    if haplotype_mode:
        command = " ".join([command, "-h"])
    command = " ".join([command, inp_fa_fNm, out_fq_one_fNm, out_fq_two_fNm,])

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )


def ssake(
    inp_fq_one_fNm,
    inp_fq_two_fNm,
    out_base_fNm,
    fragment_len,
    phred_threshold=20,  # -x
    n_consec_bases=70,  # -n
    ascii_offset=33,  # -d
    min_coverage=5,  # -w
    n_ovrlap_bases=20,  # -m
    n_reads_to_call=2,  # -o
    base_ratio=0.7,  # -r
    n_bases_to_trim=0,  # -t
    contig_size=100,  # -z
    do_track_cvrg=0,  # -c
    do_ignore_mppng=0,  # -y
    do_ignore_headr=0,  # -k
    do_break_ties=0,  # -q
    do_run_verbose=0,
):  # -v
    """
    Run a simple pipeline using SSAKE to assemble reads in a docker
    container.  SSAKE is a genomics application for de novo assembly
    of millions of very short DNA sequences. It is an easy-to-use,
    robust, reliable and tractable assembly algorithm for short
    sequence reads, such as those generated by Illumina Ltd. See:
    https://github.com/ralatsdc/SSAKE/blob/master/readme.md.

    Documentation for SSAKE version 4.0

    Usage: runSSAKE.sh read1.fq read2.fq libraryFragmentLength basename

    OPTIONS for TQSfastq
    -f   File of filenames corresponding to fasta/fastq files with reads to
         interrogate
              ! Implicitly used in this script
    -q   Phred quality score threshold (bases less than -q XX will be
         clipped, default -q 10, optional)
              ! Use -x for this script
    -n   Number of consecutive -q 10 bases (default -n 30, optional)
         TODO: Check that -n should read -e below
    -e   ASCII offset (33=standard 64=illumina, default -n 33, optional)
              ! Use -d for this script
    -v   Runs in verbose mode (-v 1 = yes, default = no, optional)

    OPTIONS for SSAKE
    -f File containing all the [paired (-p 1)] reads (required)
              ! Implicitly used in this script
              With -p 1:
                   ! Target insert size must be indicated at the end of
                     the header line (e.g. :400 for a 400bp
                     fragment/insert size)
                   ! Paired reads must be separated by ":"
                        >header:400 (or >header_barcode:400)
                        ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:GCGCACGACGCAGCACAGCAGCAGACGAC
                   Scaffolding options considered with -p 1:
                        -e Error (%) allowed on mean distance e.g. -e 0.75
                           == distance +/- 75% (default -e 0.75, optional)
                        -l Minimum number of links (read pairs) to compute
                           scaffold (default -k 5, optional)
                        -a Maximum link ratio between two best contig
                           pairs *higher values lead to least accurate
                           scaffolding* (default -a 0.3, optional)
    -g   Fasta file containing unpaired sequence reads (optional)
              ! Implicitly used in this script
    -w   Minimum depth of coverage allowed for contigs (e.g. -w 1 = process
         all reads [v3.7 behavior], required, recommended -w 5)
              ! The assembly will stop when 50+ contigs with coverage < -w have been seen.
    -s   Fasta file containing sequences to use as seeds exclusively
         (specify only if different from read set, optional)
              ! Ignored by this script
              TODO: Understand if these options are used only when -s is used
              -i Independent (de novo) assembly i.e Targets used to
                 recruit reads for de novo assembly, not guide/seed
                 reference-based assemblies (-i 1 = yes (default), 0 = no,
                 optional)
              -j Target sequence word size to hash (default -j 15)
              -u Apply read space restriction to seeds while -s option in
                 use (-u 1 = yes, default = no, optional)
    -m   Minimum number of overlapping bases with the seed/contig during
         overhang consensus build up (default -m 20)
    -o   Minimum number of reads needed to call a base during an extension
         (default -o 2)
    -r   Minimum base ratio used to accept a overhang consensus base
         (default -r 0.7)
    -t   Trim up to -t base(s) on the contig end when all possibilities have
         been exhausted for an extension (default -t 0, optional)
    -c   Track base coverage and read position for each contig (default -c 0, optional)
    -y   Ignore read mapping to consensus (-y 1 = yes, default = no, optional)
    -h   Ignore read name/header *will use less RAM if set to -h 1* (-h 1 = yes, default = no,
         optional)
              ! Use -k for this script
    -b   Base name for your output files (optional)
              ! Implicitly used in this script
    -z   Minimum contig size to track base coverage and read position
         (default -z 100, optional)
    -q   Break tie when no consensus base at position, pick random base
         (-q 1 = yes, default = no, optional)
    -p   Paired-end reads used? (-p 1 = yes, default = no, optional)
              ! Implicitly used in this script
    -v   Runs in verbose mode (-v 1 = yes, default = no, optional)
    """
    # Define image, and Docker run parameters
    image = "ralatsdio/ssake:v4.0.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define SSAKE command
    command = " ".join(
        [
            "runSSAKE.sh",
            "-x",
            str(phred_threshold),
            "-n",
            str(n_consec_bases),
            "-d",
            str(ascii_offset),
            "-w",
            str(min_coverage),
            "-m",
            str(n_ovrlap_bases),
            "-o",
            str(n_reads_to_call),
            "-r",
            str(base_ratio),
            "-t",
            str(n_bases_to_trim),
            "-z",
            str(contig_size),
            "-c",
            str(do_track_cvrg),
            "-y",
            str(do_ignore_mppng),
            "-k",
            str(do_ignore_headr),
            "-q",
            str(do_break_ties),
            "-v",
            str(do_run_verbose),
            inp_fq_one_fNm,
            inp_fq_two_fNm,
            str(fragment_len),
            out_base_fNm,
        ]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )
    return command


def spades(
    inp_fq_one_fNm,
    inp_fq_two_fNm,
    out_dir,
    *args,
    trusted_contigs_fNm=None,
    cov_cutoff="off",
    phred_offset="auto",
    **kwargs
):
    """
    Run SPAdes genome assembler v3.13.1 in a docker container.

    Usage: /home/biodocker/bin/spades.py [options] -o <output_dir>

    Basic options:
    -o <output_dir>                directory to store all the resulting files
                                       (required)
    --sc                           this flag is required for MDA (single-cell)
                                       data
    --meta                         this flag is required for metagenomic
                                       sample data
    --rna                          this flag is required for RNA-Seq
                                       data 
    --plasmid                      runs plasmidSPAdes pipeline for plasmid
                                       detection
    --iontorrent                   this flag is required for IonTorrent data
    --test                         runs SPAdes on toy dataset
    -h/--help                      prints this usage message
    -v/--version                   prints version

    Input data:
    --12 <filename>                file with interlaced forward and reverse
                                       paired-end reads
    -1 <filename>                  file with forward paired-end reads
    -2 <filename>                  file with reverse paired-end reads
    -s <filename>                  file with unpaired reads
    --merged <filename>            file with merged forward and reverse
                                       paired-end reads
    --pe<#>-12 <filename>          file with interlaced reads for paired-end
                                       library number <#> (<#> = 1,2,...,9)
    --pe<#>-1 <filename>           file with forward reads for paired-end
                                       library number <#> (<#> = 1,2,...,9)
    --pe<#>-2 <filename>           file with reverse reads for paired-end
                                       library number <#> (<#> = 1,2,...,9)
    --pe<#>-s <filename>           file with unpaired reads for paired-end
                                       library number <#> (<#> = 1,2,...,9)
    --pe<#>-m <filename>           file with merged reads for paired-end
                                       library number <#> (<#> = 1,2,...,9)
    --pe<#>-<or>                   orientation of reads for paired-end library
                                       number <#> (<#> = 1,2,...,9;
                                       <or> = fr, rf, ff)
    --s<#> <filename>              file with unpaired reads for single reads
                                       library number <#> (<#> = 1,2,...,9)
    --mp<#>-12 <filename>          file with interlaced reads for mate-pair
                                       library number <#> (<#> = 1,2,..,9)
    --mp<#>-1 <filename>           file with forward reads for mate-pair
                                       library number <#> (<#> = 1,2,..,9)
    --mp<#>-2 <filename>           file with reverse reads for mate-pair
                                       library number <#> (<#> = 1,2,..,9)
    --mp<#>-s <filename>           file with unpaired reads for mate-pair
                                       library number <#> (<#> = 1,2,..,9)
    --mp<#>-<or>                   orientation of reads for mate-pair library
                                       number <#> (<#> = 1,2,..,9;
                                       <or> = fr, rf, ff)
    --hqmp<#>-12 <filename>        file with interlaced reads for high-quality
                                       mate-pair library number <#>
                                       (<#> = 1,2,..,9)
    --hqmp<#>-1 <filename>         file with forward reads for high-quality
                                       mate-pair library number <#>
                                       (<#> = 1,2,..,9)
    --hqmp<#>-2 <filename>         file with reverse reads for high-quality
                                       mate-pair library number <#>
                                       (<#> = 1,2,..,9)
    --hqmp<#>-s <filename>         file with unpaired reads for high-quality
                                       mate-pair library number <#>
                                       (<#> = 1,2,..,9)
    --hqmp<#>-<or>                 orientation of reads for high-quality
                                       mate-pair library number <#>
                                       (<#> = 1,2,..,9; <or> = fr, rf, ff)
    --nxmate<#>-1 <filename>       file with forward reads for Lucigen NxMate
                                       library number <#> (<#> = 1,2,..,9)
    --nxmate<#>-2 <filename>       file with reverse reads for Lucigen NxMate
                                       library number <#> (<#> = 1,2,..,9)
    --sanger <filename>            file with Sanger reads
    --pacbio <filename>            file with PacBio reads
    --nanopore <filename>          file with Nanopore reads
    --tslr <filename>              file with TSLR-contigs
    --trusted-contigs <filename>   file with trusted contigs
    --untrusted-contigs <filename> file with untrusted contigs

    Pipeline options:
    --only-error-correction        runs only read error correction (without
                                       assembling)
    --only-assembler               runs only assembling (without read error
                                       correction)
    --careful                      tries to reduce number of mismatches and
                                       short indels
    --continue                     continue run from the last available
                                       check-point
    --restart-from <cp>            restart run with updated options and from
                                       the specified check-point
                                       ('ec', 'as', 'k<int>', 'mc', 'last')
    --disable-gzip-output          forces error correction not to compress the
                                       corrected reads
    --disable-rr                   disables repeat resolution stage of
                                       assembling

    Advanced options:
    --dataset <filename>           file with dataset description in YAML format
    -t/--threads <int>             number of threads [default: 16]
    -m/--memory <int>              RAM limit for SPAdes in Gb (terminates
                                       if exceeded) [default: 250]
    --tmp-dir <dirname>            directory for temporary files
                                       [default: <output_dir>/tmp]
    -k <int,int,...>               comma-separated list of k-mer sizes
                                       (must be odd and less than 128)
                                       [default: 'auto']
    --cov-cutoff <float>           coverage cutoff value (a positive float
                                       number, or 'auto', or 'off')
                                       [default: 'off']
    --phred-offset <33 or 64>      PHRED quality offset in the input reads
                                       (33 or 64) [default: auto]
    """
    # Define image, and Docker run parameters
    image = "ralatsdio/spades:v3.13.1"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define SPAdes command
    command = " ".join(["spades.py", "-1", inp_fq_one_fNm, "-2", inp_fq_two_fNm,])
    if trusted_contigs_fNm is not None:
        command = " ".join([command, "--trusted-contigs", trusted_contigs_fNm,])
    command = " ".join(
        [
            command,
            "-o",
            out_dir,
            "--cov-cutoff",
            cov_cutoff,
            "--phred-offset",
            phred_offset,
        ]
    )
    for arg in args:
        command = " ".join([command, arg,])
    for key, val in kwargs.items():
        command = " ".join([command, "--" + key.replace("_", "-"), val,])

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )
    return command


def unicycler(
    inp_fq_one_fNm, inp_fq_two_fNm, out_dir, *args, inp_fq_lng_fNm=None, **kwargs
):
    """
    Run Unicycler assembly pipeline v0.4.8 in a docker container.

    usage: unicycler [-h] [--help_all] [--version]
                     [-1 SHORT1] [-2 SHORT2] [-s UNPAIRED] [-l LONG] -o OUT
                     [--verbosity VERBOSITY] [--min_fasta_length MIN_FASTA_LENGTH]
                     [--keep KEEP] [-t THREADS] [--mode {conservative,normal,bold}]
                     [--linear_seqs LINEAR_SEQS] [--vcf]

           __
           \ \___
            \ ___\
            //
       ____//      _    _         _                     _
     //_  //\\    | |  | |       |_|                   | |
    //  \//  \\   | |  | | _ __   _   ___  _   _   ___ | |  ___  _ __
    ||  (O)  ||   | |  | || '_ \ | | / __|| | | | / __|| | / _ \| '__|
    \\    \_ //   | |__| || | | || || (__ | |_| || (__ | ||  __/| |
     \\_____//     \____/ |_| |_||_| \___| \__, | \___||_| \___||_|
                                            __/ |
                                           |___/

    Unicycler: an assembly pipeline for bacterial genomes

    Help:
      -h, --help                           Show this help message and exit
      --help_all                           Show a help message with all
                                             program options
      --version                            Show Unicycler's version number

    Input:
      -1 SHORT1, --short1 SHORT1           FASTQ file of first short reads in
                                             each pair (required)
      -2 SHORT2, --short2 SHORT2           FASTQ file of second short reads in
                                             each pair (required)
      -s UNPAIRED, --unpaired UNPAIRED     FASTQ file of unpaired short reads
                                             (optional)
      -l LONG, --long LONG                 FASTQ or FASTA file of long reads
                                             (optional)

    Output:
      -o OUT, --out OUT                    Output directory (required)
      --verbosity VERBOSITY                Level of stdout and log file
                                             information (default: 1)
                                             0 = no stdout,
                                             1 = basic progress indicators,
                                             2 = extra info,
                                             3 = debugging info
      --min_fasta_length MIN_FASTA_LENGTH  Exclude contigs from the FASTA file
                                             which are shorter than this
                                             length (default: 100)
      --keep KEEP                          Level of file retention (default: 1)
                                             0 = only keep final files:
                                                 assembly (FASTA, GFA and log),
                                             1 = also save graphs at main
                                                 checkpoints,
                                             2 = also keep SAM (enables fast
                                                 rerun in different mode),
                                             3 = keep all temp files and save
                                                 all graphs (for debugging)
      --vcf                                Produce a VCF by mapping the short
                                             reads to the final assembly
                                             (experimental, default: do not
                                              produce a vcf file)

    Other:
      -t THREADS, --threads THREADS        Number of threads used (default: 4)
      --mode {conservative,normal,bold}    Bridging mode (default: normal)
                                             conservative = smaller contigs,
                                               lowest misassembly rate
                                             normal = moderate contig size and
                                               misassembly rate
                                             bold = longest contigs, higher
                                               misassembly rate
      --linear_seqs LINEAR_SEQS            The expected number of linear
                                             (i.e. non-circular) sequences in
                                             the underlying sequence
                                             (default: 0)
    """
    # Define image, and Docker run parameters
    image = "ralatsdio/unicycler:v0.4.8"
    hosting_dir = os.path.dirname(os.path.abspath(__file__))
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define Unicycler command
    command = " ".join(["unicycler", "-1", inp_fq_one_fNm, "-2", inp_fq_two_fNm,])
    if inp_fq_lng_fNm is not None:
        command = " ".join([command, "-l", inp_fq_lng_fNm,])
    command = " ".join([command, "-o", out_dir,])
    for arg in args:
        command = " ".join([command, arg,])
    for key, val in kwargs.items():
        command = " ".join([command, "--" + key,])

    # Run the command in the Docker image
    client = docker.from_env()
    client.containers.run(
        image, command=command, volumes=volumes, working_dir=working_dir,
    )
    return command


if __name__ == "__main__":
    pass
