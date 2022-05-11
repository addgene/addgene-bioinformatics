import logging
import math
import os
import re

from Bio import SeqIO

from Bio.SeqRecord import SeqRecord
import docker

logger = logging.getLogger(__name__)


def csc(a_seq, o_seq):
    """Finds the cyclic rotation of a_seq (or an approximation of it)
    that minimizes the blockwise q-gram distance from o_seq using
    Docker image "ralatsdio/csc:v0.1.0".

    Parameters
    ----------
    a_seq : Bio.Seq.Seq
        Thought of as the reference sequence
    o_seq : Bio.Seq.Seq
        Thought of as the candidate sequence

    Returns
    -------
    tpl : tuple
        Contains the rotated reference and candidate sequences

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
    hosting_dir = os.getcwd()
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
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
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


def get_kmc_read_format(read_file_names):
    """
    Determine a valid read format implied by the read file name
    extension(s).

    Parameters
    ----------
    read_file_names : list(str)
        list of read file names

    Returns
    -------
    valid_read_format : str
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
        Read file name

    Returns
    -------
    tuple(str, boolean)
        Valid read format: "fasta" or "fastq", and flag indicating
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

    Parameters
    ----------
    read_file_names : str
        File name of reads, or file containing file names of reads
    database_file_name : str
        File name to which results are written
    k_mer_length : int
        K-mer length (k from 1 to 256; default: 25)
    signature_length : int
        Signature length (5, 6, 7, 8, 9, 10, 11); default: 9
    count_min : int
        Exclude k-mers occurring less than <value> times (default: 2)
    max_count : int
        Maximal value of a counter (default: 255)
    count_max : int
        Exclude k-mers occurring more of than <value> times (default:
            1e9)
    canonical_form : boolean
        Flag to transform k-mers into canonical form

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

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
    hosting_dir = os.getcwd()
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
        [
            command,
            read_format,
            input_str,
            database_file_name,
            working_dir,
        ]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
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
        -ci<value> - minimum value of counter to be stored in the output file
        -cx<value> - maximum value of counter to be stored in the output file

    Dump operation parameters:
        -s - sorted output

    Example:
        kmc_tools transform db \
            reduce err_kmers -cx10 \
            reduce valid_kmers -ci11 \
            histogram histo.txt \
            dump dump.txt

    Parameters
    ----------
    inp_database_file_name : str
        Path to database generated by KMC
    operation : str
        Transform operation name
    out_database_file_name : str
        Path to output
    inp_count_min : int
        Exclude k-mers occurring less than <value> times on input
        (default: 2)
    inp_count_max : int
        Exclude k-mers occurring more of than <value> times on input
        (default: 1e9)
    out_count_min : int
        Exclude k-mers occurring less than <value> times on output
        (default: 2)
    out_count_max : int
        Exclude k-mers occurring more of than <value> times on output
        (default: 1e9)
    is_sorted : boolean
        Flag to sort output (default: False)

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Check input arguments
    if operation not in ["reduce", "histogram", "dump"]:
        raise (Exception("Operation is not implemented"))

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.getcwd()
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
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


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

    Parameters
    ----------
    inp_database_file_name_a : str
        Path to database generated by KMC
    operation,
        Set operations to be performed on inputs
    inp_database_file_name_b : str
        Path to database generated by KMC
    out_database_file_name,
        Path to output k-mer database
    inp_count_min_a : int
        Exclude k-mers occurring less than <value> times on input 'a'
        (default: 2)
    inp_count_max_a : int
        Exclude k-mers occurring more of than <value> times on input
        'a' (default: 1e9)
    inp_count_min_b : int
        Exclude k-mers occurring less than <value> times on input 'b'
        (default: 2)
    inp_count_max_b : int
        Exclude k-mers occurring more of than <value> times on input
        'b' (default: 1e9)
    out_count_min : int
        Exclude k-mers occurring less than <value> times on output
        (default: 2)
    out_max_count : int
        Maximal value of a counter on output (default: 255)
    out_count_max : int
        Exclude k-mers occurring more of than <value> times on output
        (default: 1e9)
    out_calc_mode : str
        Redefine counter calculation mode for equal k-mers

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

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
    hosting_dir = os.getcwd()
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
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


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
                min, max, diff, sum, left, right (detailed description
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

    Returns
    -------
    None

    Returns
    -------
    None

    """
    raise (NotImplementedError("utility.kmc_complex() has not been implemented"))


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

    Parameters
    ----------
    inp_database_file_name : str
        Path to database generated by KMC
    inp_read_file_name : str
        Path to input set of reads
    out_read_file_name : str
        Path to output set of reads
    trim_reads : boolean
        Trim reads on first invalid k-mer instead of remove entirely
        (default: False)
    hard_mask : boolean
        Hard mask invalid k-mers in a read (default: False)
    inp_db_count_min : int
        Exclude k-mers occurring less than <value> times in database
        (default: 2)
    inp_db_count_max : int
        Exclude k-mers occurring more than <value> times in database
        (default: 1e9)
    inp_rd_count_min : int
        Remove reads containing less k-mers than value. It can be
        integer or floating number in range [0.0;1.0] (default: 2)
    inp_rd_count_max : int
        Remove reads containing more k-mers than value. It can be
        integer or floating number in range [0.0;1.0] (default: 1e9)

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Determine the input file format
    read_format = get_kmc_read_format([inp_read_file_name])

    # Define image, and Docker run parameters
    image = "ralatsdio/kmc:v3.1.1"
    hosting_dir = os.getcwd()
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
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def ssake(
    inp_fq_one_fnm,
    inp_fq_two_fnm,
    out_base_fnm,
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
    container.

    SSAKE is a genomics application for de novo assembly of millions
    of very short DNA sequences. It is an easy-to-use, robust,
    reliable and tractable assembly algorithm for short sequence
    reads, such as those generated by Illumina Ltd. See:
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

    Parameters
    ----------
    inp_fq_one_fnm : str
        File name containing read one of paired reads
    inp_fq_two_fnm : str
        File name containing read two of paired reads
    out_base_fnm : str
        Base name for your output files
    fragment_len : int
        Target insert size [bp]
    phred_threshold : int  # -x
        Phred quality score threshold (bases less than -x XX will be
        clipped) (default: 20)
    n_consec_bases : int  # -n
        Number of consecutive -x XX bases (default: 70)
    ascii_offset : int  # -d
       ASCII offset (33=standard 64=illumina) (default: 33)
    min_coverage : int  # -w
        Minimum depth of coverage allowed for contigs (default: 5)
    n_ovrlap_bases : int  # -m
        Minimum number of overlapping bases with the seed/contig
        during overhang consensus build up (default: 20)
    n_reads_to_call : int  # -o
        Minimum number of reads needed to call a base during an
        extension (default: 2)
    base_ratio : float  # -r
        Minimum base ratio used to accept a overhang consensus base
        (default: 0.7)
    n_bases_to_trim : int  # -t
        Trim up to -t base(s) on the contig end when all possibilities
        have been exhausted for an extension (default: 0)
    contig_size : int  # -z
        Minimum contig size to track base coverage and read position
        (default: 100)
    do_track_cvrg : int  # -c
        Track base coverage and read position for each contig
        (default: 0)
    do_ignore_mppng : int  # -y
        Ignore read mapping to consensus (-y 1 = yes) (default: no)
    do_ignore_headr : int  # -k
        Ignore read name/header *will use less RAM if set to -k 1* (-k
        1 = yes) (default: no)
    do_break_ties : int  # -q
        Break tie when no consensus base at position, pick random base
        (-q 1 = yes) (default: no)
    do_run_verbose : int
        Runs in verbose mode (-v 1 = yes) (default: no)

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define image, and Docker run parameters
    image = "ralatsdio/ssake:v4.0.1"
    hosting_dir = os.getcwd()
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
            inp_fq_one_fnm,
            inp_fq_two_fnm,
            str(fragment_len),
            out_base_fnm,
        ]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def spades(inp_fq_one_fnm, inp_fq_two_fnm, out_dir, *args, **kwargs):
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

    Parameters
    ----------
    inp_fq_one_fnm : str  # -1
        File with forward paired-end reads
    inp_fq_two_fnm : str  # -2
        File with reverse paired-end reads
    out_dir : str  # -o
        Directory to store all the resulting files
    args : tuple
        Values to be added as is to the command, each separated by a
        space from each other and the command
    kwargs : dict
        Key-value pairs to be added to the command, appending a double
        dash to the key, and replacing dash with underscore in the
        key, each separated by a space from each other and the
        command. Use a kwarg to provide merged reads, and trusted or
        untrusted contigs.

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define image, and Docker run parameters
    image = "ralatsdio/spades:v3.13.1"
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define SPAdes command
    command = " ".join(
        [
            "spades.py",
            "-1",
            inp_fq_one_fnm,
            "-2",
            inp_fq_two_fnm,
        ]
    )
    command = " ".join(
        [
            command,
            "-o",
            out_dir,
        ]
    )
    for arg in args:
        command = " ".join(
            [
                command,
                arg,
            ]
        )
    for key, val in kwargs.items():
        command = " ".join(
            [
                command,
                "--" + key.replace("_", "-"),
                str(val),
            ]
        )

    # Run the command in the Docker image
    client = docker.from_env()
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def apc(base_file_name, contigs_file_name):
    """
    Run apc.pl script v0.1.0 in a docker container.

    usage: /home/biodocker/bin/apc.pl [options] <fasta>

    Have you assembled [a] [p]erfect [c]ircle (e.g. plasmid, bacterial
    chromosome)? This script tests each sequence in the supplied fasta
    file for self-overlap. If overlap is found, the 5' copy of the
    overlapping sequence is trimmed, the ends joined, and the join
    moved ("permuted") to the middle of the output sequence. Join
    location is appended to the sequence header, to allow
    confirmation. If a multi-fasta file is provided, each sequence is
    tested for circularity and output separately. If reads are
    provided (e.g. PacBio corrected.fastq), BWA MEM is used to align
    reads to spliced sequence, for validation.

    The binaries 'lastdb', 'lastal', and (with -r) 'bwa' and
    'samtools' must be in your PATH.

    -b <text>   basename for output files (alphanumeric and underscore
                characters)
    -r <fastq>	high accuracy long reads to use for validation

    Parameters
    ----------
    base_file_name : str
        Basename for output files
    contigs_file_name : str
        File containing contigs to circularize

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define image, and Docker run parameters
    image = "ralatsdio/apc:v0.1.0"
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}

    # Define apc command
    command = " ".join(
        [
            "apc.pl",
            "-b",
            base_file_name,
            contigs_file_name,
        ]
    )

    # Run the command in the Docker image
    client = docker.from_env()
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def BBDuk(
    inp_fq_one_fnm,
    inp_fq_two_fnm,
    outp_fnm="output1.fastq",
    outp2_fnm="output2.fastq",
    ktrimright="t",
    k="27",
    hdist="1",
    edist="0",
    ref="adapters.fa",
    qtrim="rl",
    trimq="25",
    minlength="30",
    trimbyoverlap="t",
    minoverlap="24",
    ordered="t",
    qin="33",
    *args,
    **kwargs,
):
    """
    BBMap compares reads to the kmers in a reference dataset,
    optionally allowing an edit distance. Splits the reads into two
    outputs - those that match the reference, and those that
    don't. Can also trim (remove) the matching parts of the reads
    rather than binning the reads.  Please read
    bbmap/docs/guides/BBDukGuide.txt for more information.

    Usage:  bbduk.sh in=<input file> out=<output file> ref=<contaminant files>

    Input may be stdin or a fasta or fastq file, compressed or uncompressed.
    If you pipe via stdin/stdout, please include the file type; e.g. for gzipped
        fasta input, set in=stdin.fa.gz

    Input parameters:
        in=<file>              Main input. in=stdin.fq will pipe from stdin.
        in2=<file>             Input for 2nd read of pairs in a different file.
        ref=<file,file>        Comma-delimited list of reference files.
                                    In addition to filenames, you may also use the keywords:
                                    adapters, artifacts, phix, lambda, pjet, mtst, kapa
        literal=<seq,seq>      Comma-delimited list of literal reference sequences.
        touppercase=f          (tuc) Change all bases upper-case.
        interleaved=auto       (int) t/f overrides interleaved autodetection.
        qin=auto               Input quality offset: 33 (Sanger), 64, or auto.
        reads=-1               If positive, quit after processing X reads or pairs.
        copyundefined=f        (cu) Process non-AGCT IUPAC reference bases by making all
                                    possible unambiguous copies.  Intended for short motifs
                                    or adapter barcodes, as time/memory use is exponential.
        samplerate=1           Set lower to only process a fraction of input reads.
        samref=<file>          Optional reference fasta for processing sam files.

    Output parameters:
        out=<file>             (outnonmatch) Write reads here that do not contain
                                    kmers matching the database.  'out=stdout.fq' will pipe
                                    to standard out.
        out2=<file>            (outnonmatch2) Use this to write 2nd read of pairs to a
                                   different file.
        outm=<file>            (outmatch) Write reads here that fail filters.  In default
                                   kfilter mode, this means any read with a matching kmer.
                                   In any mode, it also includes reads that fail filters such
                                   as minlength, mingc, maxgc, entropy, etc.  In other words,
                                   it includes all reads that do not go to 'out'.
        outm2=<file>           (outmatch2) Use this to write 2nd read of pairs to a
                                   different file.
        outs=<file>            (outsingle) Use this to write singleton reads whose mate
                                   was trimmed shorter than minlen.
        stats=<file>           Write statistics about which contaminants were detected.
        refstats=<file>        Write statistics on a per-reference-file basis.
        rpkm=<file>            Write RPKM for each reference sequence (for RNA-seq).
        dump=<file>            Dump kmer tables to a file, in fasta format.
        duk=<file>             Write statistics in duk's format. *DEPRECATED*
        nzo=t                  Only write statistics about ref sequences with nonzero hits.
        overwrite=t            (ow) Grant permission to overwrite files.
        showspeed=t           (ss) 'f' suppresses display of processing speed.
        ziplevel=2             (zl) Compression level; 1 (min) through 9 (max).
        fastawrap=70           Length of lines in fasta output.
        qout=auto              Output quality offset: 33 (Sanger), 64, or auto.
        statscolumns=3         (cols) Number of columns for stats output, 3 or 5.
                                   5 includes base counts.
        rename=f               Rename reads to indicate which sequences they matched.
        refnames=f             Use names of reference files rather than scaffold IDs.
        trd=f                  Truncate read and ref names at the first whitespace.
        ordered=f              Set to true to output reads in same order as input.
        maxbasesout=-1         If positive, quit after writing approximately this many
                                   bases to out (outu/outnonmatch).
        maxbasesoutm=-1        If positive, quit after writing approximately this many
                                   bases to outm (outmatch).
        json=f                 Print to screen in json format.

    Histogram output parameters:
        bhist=<file>          Base composition histogram by position.
        qhist=<file>          Quality histogram by position.
        qchist=<file>         Count of bases with each quality value.
        aqhist=<file>         Histogram of average read quality.
        bqhist=<file>         Quality histogram designed for box plots.
        lhist=<file>          Read length histogram.
        phist=<file>          Polymer length histogram.
        gchist=<file>         Read GC content histogram.
        ihist=<file>          Insert size histogram, for paired reads in mapped sam.
        gcbins=100            Number gchist bins.  Set to 'auto' to use read length.
        maxhistlen=6000       Set an upper bound for histogram lengths; higher uses
                                more memory.  The default is 6000 for some histograms
                                and 80000 for others.

    Histograms for mapped sam/bam files only:
        histbefore=t         Calculate histograms from reads before processing.
        ehist=<file>         Errors-per-read histogram.
        qahist=<file>        Quality accuracy histogram of error rates versus quality
        score.
        indelhist=<file>     Indel length histogram.
        mhist=<file>         Histogram of match, sub, del, and ins rates by position.
        idhist=<file>        Histogram of read count versus percent identity.
        idbins=100           Number idhist bins.  Set to 'auto' to use read length.
        varfile=<file>       Ignore substitution errors listed in this file when
                                 calculating error rates.  Can be generated with
                                 CallVariants.
        vcf=<file>           Ignore substitution errors listed in this VCF file
                                 when calculating error rates.
        ignorevcfindels=t    Also ignore indels listed in the VCF.

    Processing parameters:
        k=27                Kmer length used for finding contaminants.  Contaminants
                                shorter than k will not be found.  k must be at least 1.
        rcomp=t             Look for reverse-complements of kmers in addition to
                                forward kmers.
        maskmiddle=t        (mm) Treat the middle base of a kmer as a wildcard, to
                                increase sensitivity in the presence of errors.
        minkmerhits=1       (mkh) Reads need at least this many matching kmers
                                to be considered as matching the reference.
        minkmerfraction=0.0 (mkf) A reads needs at least this fraction of its total
                                kmers to hit a ref, in order to be considered a match.
                                If this and minkmerhits are set, the greater is used.
        mincovfraction=0.0  (mcf) A reads needs at least this fraction of its total
                                bases to be covered by ref kmers to be considered a match.
                                If specified, mcf overrides mkh and mkf.
        hammingdistance=0   (hdist) Maximum Hamming distance for ref kmers (subs only).
                                Memory use is proportional to (3*K)^hdist.
        qhdist=0            Hamming distance for query kmers; impacts speed, not memory.
        editdistance=0      (edist) Maximum edit distance from ref kmers (subs
                                                                          and indels).  Memory use is proportional to (8*K)^edist.
        hammingdistance2=0  (hdist2) Sets hdist for short kmers, when using mink.
        qhdist2=0           Sets qhdist for short kmers, when using mink.
        editdistance2=0     (edist2) Sets edist for short kmers, when using mink.
        forbidn=f           (fn) Forbids matching of read kmers containing N.
                                By default, these will match a reference 'A' if
                                hdist>0 or edist>0, to increase sensitivity.
        removeifeitherbad=t (rieb) Paired reads get sent to 'outmatch' if either is
                                match (or either is trimmed shorter than minlen).
                                Set to false to require both.
        trimfailures=f      Instead of discarding failed reads, trim them to 1bp.
                                This makes the statistics a bit odd.
        findbestmatch=f     (fbm) If multiple matches, associate read with sequence
                                sharing most kmers.  Reduces speed.
        skipr1=f            Don't do kmer-based operations on read 1.
        skipr2=f            Don't do kmer-based operations on read 2.
        ecco=f              For overlapping paired reads only.  Performs error-
                                correction with BBMerge prior to kmer operations.
                                recalibrate=f       (recal) Recalibrate quality scores.  Requires calibration
                                matrices generated by CalcTrueQuality.
        sam=<file,file>     If recalibration is desired, and matrices have not already
                                been generated, BBDuk will create them from the sam file.
        amino=f             Run in amino acid mode.  Some features have not been
                                tested, but kmer-matching works fine.  Maximum k is 12.

    Speed and Memory parameters:
        threads=auto        (t) Set number of threads to use; default is number of
                                logical processors.
        prealloc=f          Preallocate memory in table.  Allows faster table loading
                                and more efficient memory usage, for a large reference.
                                monitor=f           Kill this process if it crashes.  monitor=600,0.01 would
                                kill after 600 seconds under 1% usage.
        minrskip=1          (mns) Force minimal skip interval when indexing reference
                                kmers.  1 means use all, 2 means use every other kmer, etc.
        maxrskip=1          (mxs) Restrict maximal skip interval when indexing
                                reference kmers. Normally all are used for scaffolds<100kb,
                                but with longer scaffolds, up to maxrskip-1 are skipped.
        rskip=              Set both minrskip and maxrskip to the same value.
                                If not set, rskip will vary based on sequence length.
        qskip=1             Skip query kmers to increase speed.  1 means use all.
        speed=0             Ignore this fraction of kmer space (0-15 out of 16) in both
                                reads and reference.  Increases speed and reduces memory.
        Note: Do not use more than one of 'speed', 'qskip', and 'rskip'.

    Trimming/Filtering/Masking parameters:
        Note - if ktrim, kmask, and ksplit are unset, the default behavior is kfilter.
        All kmer processing modes are mutually exclusive.
        Reads only get sent to 'outm' purely based on kmer matches in kfilter mode.

        ktrim=f             Trim reads to remove bases matching reference kmers.
                                Values:
                                f (don't trim),
                                r (trim to the right),
                                l (trim to the left)
        kmask=              Replace bases matching ref kmers with another symbol.
                                Allows any non-whitespace character, and processes short
                                kmers on both ends if mink is set.  'kmask=lc' will
                                convert masked bases to lowercase.
        maskfullycovered=f  (mfc) Only mask bases that are fully covered by kmers.
        ksplit=f            For single-ended reads only.  Reads will be split into
                                pairs around the kmer.  If the kmer is at the end of the
                                read, it will be trimmed instead.  Singletons will go to
                                out, and pairs will go to outm.  Do not use ksplit with
                                other operations such as quality-trimming or filtering.
                                mink=0              Look for shorter kmers at read tips down to this length,
                                when k-trimming or masking.  0 means disabled.  Enabling
                                this will disable maskmiddle.
        qtrim=f             Trim read ends to remove bases with quality below trimq.
                                Performed AFTER looking for kmers.  Values:
                                    rl (trim both ends),
                                f (neither end),
                                r (right end only),
                                l (left end only),
                                w (sliding window).
        trimq=6             Regions with average quality BELOW this will be trimmed,
                                if qtrim is set to something other than f.  Can be a
                                floating-point number like 7.3.
        trimclip=f          Trim soft-clipped bases from sam files.
                                minlength=10        (ml) Reads shorter than this after trimming will be
                                discarded.  Pairs will be discarded if both are shorter.
                                mlf=0               (minlengthfraction) Reads shorter than this fraction of
                                original length after trimming will be discarded.
                                maxlength=          Reads longer than this after trimming will be discarded.
        minavgquality=0     (maq) Reads with average quality (after trimming) below
                                this will be discarded.
        maqb=0              If positive, calculate maq from this many initial bases.
        minbasequality=0    (mbq) Reads with any base below this quality (after
                                trimming) will be discarded.
        maxns=-1            If non-negative, reads with more Ns than this
                                (after trimming) will be discarded.
        mcb=0               (minconsecutivebases) Discard reads without at least
                                this many consecutive called bases.
        ottm=f              (outputtrimmedtomatch) Output reads trimmed to shorter
                                than minlength to outm rather than discarding.
        tp=0                (trimpad) Trim this much extra around matching kmers.
        tbo=f               (trimbyoverlap) Trim adapters based on where paired
                                reads overlap.
        strictoverlap=t     Adjust sensitivity for trimbyoverlap mode.
        minoverlap=14       Require this many bases of overlap for detection.
        mininsert=40        Require insert size of at least this for overlap.
                                Should be reduced to 16 for small RNA sequencing.
        tpe=f               (trimpairsevenly) When kmer right-trimming, trim both
                                reads to the minimum length of either.
        forcetrimleft=0     (ftl) If positive, trim bases to the left of this position
                                (exclusive, 0-based).
        forcetrimright=0    (ftr) If positive, trim bases to the right of this position
                                (exclusive, 0-based).
        forcetrimright2=0   (ftr2) If positive, trim this many bases on the right end.
        forcetrimmod=0      (ftm) If positive, right-trim length to be equal to zero,
                                modulo this number.
        restrictleft=0      If positive, only look for kmer matches in the
                                leftmost X bases.
        restrictright=0     If positive, only look for kmer matches in the
                                rightmost X bases.
        mingc=0             Discard reads with GC content below this.
        maxgc=1             Discard reads with GC content above this.
        gcpairs=t           Use average GC of paired reads.
                                Also affects gchist.
        tossjunk=f          Discard reads with invalid characters as bases.
        swift=f             Trim Swift sequences: Trailing C/T/N R1, leading G/A/N R2.

                                                                                       Header-parsing parameters - these require Illumina headers:
        chastityfilter=f    (cf) Discard reads with id containing ' 1:Y:' or ' 2:Y:'.
        barcodefilter=f     Remove reads with unexpected barcodes if barcodes is set,
                                or barcodes containing 'N' otherwise.  A barcode must be
                                the last part of the read header.  Values:
                                t:     Remove reads with bad barcodes.
                                f:     Ignore barcodes.
        crash: Crash upon encountering bad barcodes.
        barcodes=           Comma-delimited list of barcodes or files of barcodes.
        xmin=-1             If positive, discard reads with a lesser X coordinate.
        ymin=-1             If positive, discard reads with a lesser Y coordinate.
        xmax=-1             If positive, discard reads with a greater X coordinate.
        ymax=-1             If positive, discard reads with a greater Y coordinate.

    Polymer trimming:
        trimpolya=0         If greater than 0, trim poly-A or poly-T tails of
                                at least this length on either end of reads.
        trimpolygleft=0     If greater than 0, trim poly-G prefixes of at least this
                                length on the left end of reads.  Does not trim poly-C.
        trimpolygright=0    If greater than 0, trim poly-G tails of at least this
                                length on the right end of reads.  Does not trim poly-C.
        trimpolyg=0         This sets both left and right at once.
        filterpolyg=0       If greater than 0, remove reads with a poly-G prefix of
                                at least this length (on the left).
        Note: there are also equivalent poly-C flags.

    Polymer tracking:
        pratio=base,base    'pratio=G,C' will print the ratio of G to C polymers.
        plen=20             Length of homopolymers to count.

                                                          Entropy/Complexity parameters:
        entropy=-1          Set between 0 and 1 to filter reads with entropy below
                                that value.  Higher is more stringent.
        entropywindow=50    Calculate entropy using a sliding window of this length.
        entropyk=5          Calculate entropy using kmers of this length.
        minbasefrequency=0  Discard reads with a minimum base frequency below this.
        entropytrim=f       Values:
                                f:  (false) Do not entropy-trim.
                                r:  (right) Trim low entropy on the right end only.
                                l:  (left) Trim low entropy on the left end only.
                                rl: (both) Trim low entropy on both ends.
        entropymask=f       Values:
                                f:  (filter) Discard low-entropy sequences.
                                t:  (true) Mask low-entropy parts of sequences with N.
                                lc: Change low-entropy parts of sequences to lowercase.
        entropymark=f       Mark each base with its entropy value.  This is on a scale
                                of 0-41 and is reported as quality scores, so the output
                                should be fastq or fasta+qual.
        NOTE: If set, entropytrim overrides entropymask.

    Cardinality estimation:
        cardinality=f       (loglog) Count unique kmers using the LogLog algorithm.
        cardinalityout=f    (loglogout) Count unique kmers in output reads.
        loglogk=31          Use this kmer length for counting.
        loglogbuckets=1999  Use this many buckets for counting.

    Java Parameters:

        -Xmx                This will set Java's memory usage, overriding autodetection.
                                -Xmx20g will
                                specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                                The max is typically 85% of physical memory.
        -eoom               This flag will cause the process to exit if an
                                out-of-memory exception occurs.  Requires Java 8u92+.
        -da                 Disable assertions.


    Parameters
    ----------
    inp_fq_one_fnm : str
        Main input
    inp_fq_two_fnm : str
        Input for 2nd read of pairs in a different file
    outp_fnm : str
        Write reads here that do not contain kmers matching the
        database (default: 'output1.fastq')
    outp2_fnm : str
        Use this to write 2nd read of pairs to a different file
        (default: 'output2.fastq')
    ktrimright : str
        Flag to trim to the right (default: 't')
    k : int
        Kmer length used for finding contaminants (default: 27)
    hdist : int
        Maximum Hamming distance for ref kmers (subs only) (default:
        1)
    edist : int
        Maximum edit distance from ref kmers (subs and indels)
        (default: 0)
    ref : str
        Name of reference file (default: 'adapters.fa')
    qtrim : str
        Trim neither, right, left or sliding window read ends to
        remove bases with quality below trimq (default: 'rl')
    trimq : int
        Regions with average quality BELOW this will be trimmed, if
        qtrim is set to something other than f (default: 25)
    minlength : int
        Reads shorter than this after trimming will be discarded.
        Pairs will be discarded if both are shorter. (default: 30)
    trimbyoverlap : str
        Trim adapters based on where paired reads overlap (default:
        't')
    minoverlap : int
        Require this many bases of overlap for detection (default: 24)
    ordered : str
        Set to true to output reads in same order as input (default:
        't')
    qin : int
        Input quality offset: 33 (Sanger), 64, or auto (default: 33)
    args : tuple
        Values to be added as is to the command, each separated by a
        space from each other and the command
    kwargs : dict
        Key-value pairs to be added to the command, appending a double
        dash to the key, and replacing dash with underscore in the
        key, each separated by a space from each other and the
        command. Use a kwarg to provide trusted or untrusted contigs

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define BBDuk command
    command = " ".join(
        [
            "bbduk.sh",
            f"in={inp_fq_one_fnm}",
            f"in2={inp_fq_two_fnm}",
            f"out={outp_fnm}",
            f"out2={outp2_fnm}",
            f"ktrimright={ktrimright}",
            f"k={k}",
            f"hdist={hdist}",
            f"edist={edist}",
            f"ref={ref}",
            f"qtrim={qtrim}",
            f"trimq={trimq}",
            f"minlength={minlength}",
            f"trimbyoverlap={trimbyoverlap}",
            f"minoverlap={minoverlap}",
            f"ordered={ordered}",
            f"qin={qin}",
        ]
    )
    for arg in args:
        command = " ".join([command, arg])
    for key, val in kwargs.items():
        command = " ".join(
            [
                command,
                "--" + key.replace("_", "-"),
                str(val),
            ]
        )

    # Run the command in the Docker image
    image = "ralatsdio/bbtools:v38.90"
    client = docker.from_env()
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def BBNorm(
    inp_fq_one_fnm,
    inp_fq_two_fnm,
    outp_fnm="output1.fastq",
    outp2_fnm="output2.fastq",
    target="100",
    mindepth="6",
    threads="8",
    qin="33",
    *args,
    **kwargs,
):
    """
    BBNorm normalizes read depth based on kmer counts.  Can also
    error-correct, bin reads by kmer depth, and generate a kmer depth
    histogram.  However, Tadpole has superior error-correction to
    BBNorm.  Please read bbmap/docs/guides/BBNormGuide.txt for more
    information.

    Usage:     bbnorm.sh in=<input> out=<reads to keep> outt=<reads to toss> hist=<histogram output>

    Input parameters:
        in=null                Primary input.  Use in2 for paired reads in a second file
        in2=null               Second input file for paired reads in two files
        extra=null             Additional files to use for input (generating hash table) but not for output
        fastareadlen=2^31      Break up FASTA reads longer than this.  Can be useful when processing scaffolded genomes
        tablereads=-1          Use at most this many reads when building the hashtable (-1 means all)
        kmersample=1           Process every nth kmer, and skip the rest
        readsample=1           Process every nth read, and skip the rest
        interleaved=auto       May be set to true or false to force the input read file to override autodetection of the input file as paired interleaved.
        qin=auto               ASCII offset for input quality.  May be 33 (Sanger), 64 (Illumina), or auto.

    Output parameters:
        out=<file>             File for normalized or corrected reads.  Use out2 for paired reads in a second file
        outt=<file>            (outtoss) File for reads that were excluded from primary output
        reads=-1               Only process this number of reads, then quit (-1 means all)
        sampleoutput=t         Use sampling on output as well as input (not used if sample rates are 1)
        keepall=f              Set to true to keep all reads (e.g. if you just want error correction).
        zerobin=f              Set to true if you want kmers with a count of 0 to go in the 0 bin instead of the 1 bin in histograms.
                                   Default is false, to prevent confusion about how there can be 0-count kmers.
                                   The reason is that based on the 'minq' and 'minprob' settings, some kmers may be excluded from the bloom filter.
        tmpdir=                This will specify a directory for temp files (only needed for multipass runs).  If null, they will be written to the output directory.
        usetempdir=t           Allows enabling/disabling of temporary directory; if disabled, temp files will be written to the output directory.
        qout=auto              ASCII offset for output quality.  May be 33 (Sanger), 64 (Illumina), or auto (same as input).
        rename=f               Rename reads based on their kmer depth.

    Hashing parameters:
        k=31                   Kmer length (values under 32 are most efficient, but arbitrarily high values are supported)
        bits=32                Bits per cell in bloom filter; must be 2, 4, 8, 16, or 32.  Maximum kmer depth recorded is 2^cbits.  Automatically reduced to 16 in 2-pass. Large values decrease accuracy for a fixed amount of memory, so use the lowest number you can that will still capture highest-depth kmers.
        hashes=3               Number of times each kmer is hashed and stored.  Higher is slower.
                                   Higher is MORE accurate if there is enough memory, and LESS accurate if there is not enough memory.
        prefilter=f            True is slower, but generally more accurate; filters out low-depth kmers from the main hashtable.  The prefilter is more memory-efficient because it uses 2-bit cells.
        prehashes=2            Number of hashes for prefilter.
        prefilterbits=2        (pbits) Bits per cell in prefilter.
        prefiltersize=0.35     Fraction of memory to allocate to prefilter.
        buildpasses=1          More passes can sometimes increase accuracy by iteratively removing low-depth kmers
        minq=6                 Ignore kmers containing bases with quality below this
        minprob=0.5            Ignore kmers with overall probability of correctness below this
        threads=auto           (t) Spawn exactly X hashing threads (default is number of logical processors).  Total active threads may exceed X due to I/O threads.
        rdk=t                  (removeduplicatekmers) When true, a kmer's count will only be incremented once per read pair, even if that kmer occurs more than once.

    Normalization parameters:
        fixspikes=f            (fs) Do a slower, high-precision bloom filter lookup of kmers that appear to have an abnormally high depth due to collisions.
        target=100             (tgt) Target normalization depth.  NOTE:  All depth parameters control kmer depth, not read depth.
                                   For kmer depth Dk, read depth Dr, read length R, and kmer size K:  Dr=Dk*(R/(R-K+1))
        maxdepth=-1            (max) Reads will not be downsampled when below this depth, even if they are above the target depth.
        mindepth=5             (min) Kmers with depth below this number will not be included when calculating the depth of a read.
        minkmers=15            (mgkpr) Reads must have at least this many kmers over min depth to be retained.  Aka 'mingoodkmersperread'.
        percentile=54.0        (dp) Read depth is by default inferred from the 54th percentile of kmer depth, but this may be changed to any number 1-100.
        uselowerdepth=t        (uld) For pairs, use the depth of the lower read as the depth proxy.
        deterministic=t        (dr) Generate random numbers deterministically to ensure identical output between multiple runs.  May decrease speed with a huge number of threads.
        passes=2               (p) 1 pass is the basic mode.  2 passes (default) allows greater accuracy, error detection, better control of output depth.

    Error detection parameters:
        hdp=90.0                (highdepthpercentile) Position in sorted kmer depth array used as proxy of a read's high kmer depth.
        ldp=25.0                (lowdepthpercentile) Position in sorted kmer depth array used as proxy of a read's low kmer depth.
        tossbadreads=f          (tbr) Throw away reads detected as containing errors.
        requirebothbad=f        (rbb) Only toss bad pairs if both reads are bad.
        errordetectratio=125    (edr) Reads with a ratio of at least this much between their high and low depth kmers will be classified as error reads.
        highthresh=12           (ht) Threshold for high kmer.  A high kmer at this or above are considered non-error.
        lowthresh=3             (lt) Threshold for low kmer.  Kmers at this and below are always considered errors.

    Error correction parameters:
        ecc=f                   Set to true to correct errors.  NOTE: Tadpole is now preferred for ecc as it does a better job.
        ecclimit=3              Correct up to this many errors per read.  If more are detected, the read will remain unchanged.
        errorcorrectratio=140   (ecr) Adjacent kmers with a depth ratio of at least this much between will be classified as an error.
        echighthresh=22         (echt) Threshold for high kmer.  A kmer at this or above may be considered non-error.
        eclowthresh=2           (eclt) Threshold for low kmer.  Kmers at this and below are considered errors.
        eccmaxqual=127          Do not correct bases with quality above this value.
        aec=f                   (aggressiveErrorCorrection) Sets more aggressive values of ecr=100, ecclimit=7, echt=16, eclt=3.
        cec=f                   (conservativeErrorCorrection) Sets more conservative values of ecr=180, ecclimit=2, echt=30, eclt=1, sl=4, pl=4.
        meo=f                   (markErrorsOnly) Marks errors by reducing quality value of suspected errors; does not correct anything.
        mue=t                   (markUncorrectableErrors) Marks errors only on uncorrectable reads; requires 'ecc=t'.
        overlap=f               (ecco) Error correct by read overlap.

    Depth binning parameters:
        lowbindepth=10          (lbd) Cutoff for low depth bin.
        highbindepth=80         (hbd) Cutoff for high depth bin.
        outlow=<file>           Pairs in which both reads have a median below lbd go into this file.
        outhigh=<file>          Pairs in which both reads have a median above hbd go into this file.
        outmid=<file>           All other pairs go into this file.

    Histogram parameters:
        hist=<file>             Specify a file to write the input kmer depth histogram.
        histout=<file>          Specify a file to write the output kmer depth histogram.
        histcol=3               (histogramcolumns) Number of histogram columns, 2 or 3.
        pzc=f                   (printzerocoverage) Print lines in the histogram with zero coverage.
        histlen=1048576         Max kmer depth displayed in histogram.  Also affects statistics displayed, but does not affect normalization.

    Peak calling parameters:
        peaks=<file>            Write the peaks to this file.  Default is stdout.
        minHeight=2             (h) Ignore peaks shorter than this.
        minVolume=5             (v) Ignore peaks with less area than this.
        minWidth=3              (w) Ignore peaks narrower than this.
        minPeak=2               (minp) Ignore peaks with an X-value below this.
        maxPeak=BIG             (maxp) Ignore peaks with an X-value above this.
        maxPeakCount=8          (maxpc) Print up to this many peaks (prioritizing height).

    Java Parameters:
        -Xmx                    This will set Java's memory usage, overriding autodetection.
                                    -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                                    The max is typically 85% of physical memory.
        -eoom                   This flag will cause the process to exit if an
                                    out-of-memory exception occurs.  Requires Java 8u92+.
        -da                     Disable assertions.

    Parameters
    ----------
    inp_fq_one_fnm : str
        Primary input
    inp_fq_two_fnm : str
        Second input file for paired reads in two files
    outp_fnm : str
        File for normalized or corrected reads (default:
        'output1.fastq')
    outp2_fnm : str
        Second output file for paired reads in two files (default:
        'output2.fastq')
    target : int
        Target normalization depth.  All depth parameters control kmer
        depth, not read depth. (default: 100)
    mindepth : int
        Kmers with depth below this number will not be included when
        calculating the depth of a read (default: 6)
    threads : int
        Spawn exactly X hashing threads (default is number of logical
        processors) (default: 8)
    qin : int
        ASCII offset for input quality.  May be 33 (Sanger), 64
        (Illumina), or auto. (default: 33)
    args : tuple
        Values to be added as is to the command, each separated by a
        space from each other and the command
    kwargs : dict
        Key-value pairs to be added to the command, appending a double
        dash to the key, and replacing dash with underscore in the
        key, each separated by a space from each other and the
        command. Use a kwarg to provide trusted or untrusted contigs

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define BBNorm command
    command = " ".join(
        [
            "bbnorm.sh",
            f"in={inp_fq_one_fnm}",
            f"in2={inp_fq_two_fnm}",
            f"out={outp_fnm}",
            f"out2={outp2_fnm}",
            f"target={target}",
            f"mindepth={mindepth}",
            f"threads={threads}",
            f"qin={qin}",
        ]
    )
    for arg in args:
        command = " ".join([command, arg])
    for key, val in kwargs.items():
        command = " ".join(
            [
                command,
                "--" + key.replace("_", "-"),
                str(val),
            ]
        )

    # Run the command in the Docker image
    image = "ralatsdio/bbtools:v38.90"
    client = docker.from_env()
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def BBMerge(
    inp_fq_one_fnm,
    inp_fq_two_fnm,
    outp_fnm="merged.fastq",
    outpu1_fnm="unmerged1.fastq",
    outpu2_fnm="unmerged2.fastq",
    strictness="vloose=t",
    qin="33",
    *args,
    **kwargs,
):
    """
    BBmerge merges paired reads into single reads by overlap
    detection.  With sufficient coverage, can merge nonoverlapping
    reads by kmer extension.  Kmer modes (Tadpole or Bloom Filter)
    require much more memory, and should be used with the
    bbmerge-auto.sh script rather than bbmerge.sh.  Please read
    bbmap/docs/guides/BBMergeGuide.txt for more information.

    Usage for interleaved files:    bbmerge.sh in=<reads> out=<merged reads> outu=<unmerged reads>
    Usage for paired files:         bbmerge.sh in1=<read1> in2=<read2> out=<merged reads> outu1=<unmerged1> outu2=<unmerged2>

    Input may be stdin or a file, fasta or fastq, raw or gzipped.

    Input parameters:
        in=null                Primary input. 'in2' will specify a second file.
        interleaved=auto       May be set to true or false to override autodetection of
                               whether the input file as interleaved.
        reads=-1               Quit after this many read pairs (-1 means all).

    Output parameters:
        out=<file>             File for merged reads. 'out2' will specify a second file.
        outu=<file>            File for unmerged reads. 'outu2' will specify a second file.
        outinsert=<file>       (outi) File to write read names and insert sizes.
        outadapter=<file>      (outa) File to write consensus adapter sequences.
        outc=<file>            File to write input read kmer cardinality estimate.
        ihist=<file>           (hist) Insert length histogram output file.
        nzo=t                  Only print histogram bins with nonzero values.
        showhiststats=t        Print extra header lines with statistical information.
        ziplevel=2             Set to 1 (lowest) through 9 (max) to change compression
                                   level; lower compression is faster.
        ordered=f              Output reads in same order as input.
        mix=f                  Output both the merged (or mergable) and unmerged reads
                                   in the same file (out=).  Useful for ecco mode.

    Trimming/Filtering parameters:
        qtrim=f                Trim read ends to remove bases with quality below minq.
                                   Trims BEFORE merging.
                                   Values: t (trim both ends),
                                           f (neither end),
                                           r (right end only),
                                           l (left end only).
        qtrim2=f               May be specified instead of qtrim to perform trimming
                                   only if merging is unsuccessful, then retry merging.
        trimq=10               Trim quality threshold.  This may be a comma-delimited
                                   list (ascending) to try multiple values.
        minlength=1            (ml) Reads shorter than this after trimming, but before
                                   merging, will be discarded. Pairs will be discarded only
                                   if both are shorter.
        maxlength=-1           Reads with longer insert sizes will be discarded.
        tbo=f                  (trimbyoverlap) Trim overlapping reads to remove
                                   rightmost (3') non-overlapping portion, instead of joining.
        minavgquality=0        (maq) Reads with average quality below this, after
                                   trimming, will not be attempted to be merged.
        maxexpectederrors=0    (mee) If positive, reads with more combined expected
                                   errors than this will not be attempted to be merged.
        forcetrimleft=0        (ftl) If nonzero, trim left bases of the read to
                                   this position (exclusive, 0-based).
        forcetrimright=0       (ftr) If nonzero, trim right bases of the read
                                   after this position (exclusive, 0-based).
        forcetrimright2=0      (ftr2) If positive, trim this many bases on the right end.
        forcetrimmod=5         (ftm) If positive, trim length to be equal to
                                   zero modulo this number.
        ooi=f                  Output only incorrectly merged reads, for testing.
        trimpolya=t            Trim trailing poly-A tail from adapter output.  Only
                                   affects outadapter.  This also trims poly-A followed
                                   by poly-G, which occurs on NextSeq.

    Processing Parameters:
        usejni=f               (jni) Do overlapping in C code, which is faster.  Requires
                                   compiling the C code; details are in /jni/README.txt.
                                   However, the jni path is currently disabled.
        merge=t                Create merged reads.  If set to false, you can still
                                   generate an insert histogram.
        ecco=f                 Error-correct the overlapping part, but don't merge.
        trimnonoverlapping=f   (tno) Trim all non-overlapping portions, leaving only
                                   consensus sequence.  By default, only sequence to the
                                   right of the overlap (adapter sequence) is trimmed.
        useoverlap=t           Attempt find the insert size using read overlap.
        mininsert=35           Minimum insert size to merge reads.
        mininsert0=35          Insert sizes less than this will not be considered.
                                   Must be less than or equal to mininsert.
        minoverlap=12          Minimum number of overlapping bases to allow merging.
        minoverlap0=8          Overlaps shorter than this will not be considered.
                                   Must be less than or equal to minoverlap.
        minq=9                 Ignore bases with quality below this.
        maxq=41                Cap output quality scores at this.
        entropy=t              Increase the minimum overlap requirement for low-
                                   complexity reads.
        efilter=6              Ban overlaps with over this many times the expected
                                   number of errors.  Lower is more strict. -1 disables.
        pfilter=0.00004        Ban improbable overlaps.  Higher is more strict. 0 will
                                   disable the filter; 1 will allow only perfect overlaps.
        kfilter=0              Ban overlaps that create kmers with count below
                                   this value (0 disables).  If this is used minprob should
                                   probably be set to 0.  Requires good coverage.
        ouq=f                  Calculate best overlap using quality values.
        owq=t                  Calculate best overlap without using quality values.
        usequality=t           If disabled, quality values are completely ignored,
                                   both for overlap detection and filtering.  May be useful
                                   for data with inaccurate quality values.
        iupacton=f             (itn) Change ambiguous IUPAC symbols to N.
        adapter=               Specify the adapter sequences used for these reads, if
                                   known; this can be a fasta file or a literal sequence.
                                   Read 1 and 2 can have adapters specified independently
                                   with the adapter1 and adapter2 flags.  adapter=default
                                   will use a list of common adapter sequences.

    Ratio Mode:
        ratiomode=t            Score overlaps based on the ratio of matching to
                                   mismatching bases.
        maxratio=0.09          Max error rate; higher increases merge rate.
        ratiomargin=5.5        Lower increases merge rate; min is 1.
        ratiooffset=0.55       Lower increases merge rate; min is 0.
        maxmismatches=20       Maximum mismatches allowed in overlapping region.
        ratiominoverlapreduction=3  This is the difference between minoverlap in
                                    flat mode and minoverlap in ratio mode; generally,
                                    minoverlap should be lower in ratio mode.
        minsecondratio=0.1     Cutoff for second-best overlap ratio.
        forcemerge=f           Disable all filters and just merge everything
                                 (not recommended).

    Flat Mode:
        flatmode=f             Score overlaps based on the total number of mismatching
                                   bases only.
        margin=2               The best overlap must have at least 'margin' fewer
                                   mismatches than the second best.
        mismatches=3           Do not allow more than this many mismatches.
        requireratiomatch=f    (rrm) Require the answer from flat mode and ratio mode
                                   to agree, reducing false positives if both are enabled.
        trimonfailure=t        (tof) If detecting insert size by overlap fails,
                                   the reads will be trimmed and this will be re-attempted.


    *** Ratio Mode and Flat Mode may be used alone or simultaneously. ***
    *** Ratio Mode is usually more accurate and is the default mode. ***


    Strictness (these are mutually exclusive macros that set other parameters):
        strict=f               Decrease false positive rate and merging rate.
        verystrict=f           (vstrict) Greatly decrease FP and merging rate.
        ultrastrict=f          (ustrict) Decrease FP and merging rate even more.
        maxstrict=f            (xstrict) Maximally decrease FP and merging rate.
        loose=f                Increase false positive rate and merging rate.
        veryloose=f            (vloose) Greatly increase FP and merging rate.
        ultraloose=f           (uloose) Increase FP and merging rate even more.
        maxloose=f             (xloose) Maximally decrease FP and merging rate.
        fast=f                 Fastest possible mode; less accurate.

    Tadpole Parameters (for read extension and error-correction):
    *Note: These require more memory and should be run with bbmerge-auto.sh.*
        k=31                   Kmer length.  31 (or less) is fastest and uses the least
                                   memory, but higher values may be more accurate.
                                   60 tends to work well for 150bp reads.
        extend=0               Extend reads to the right this much before merging.
                                   Requires sufficient (>5x) kmer coverage.
        extend2=0              Extend reads this much only after a failed merge attempt,
                                   or in rem/rsem mode.
        iterations=1           (ei) Iteratively attempt to extend by extend2 distance
                                   and merge up to this many times.
        rem=f                  (requireextensionmatch) Do not merge if the predicted
                                   insert size differs before and after extension.
                                   However, if only the extended reads overlap, then that
                                   insert will be used.  Requires setting extend2.
        rsem=f                 (requirestrictextensionmatch) Similar to rem but stricter.
                                   Reads will only merge if the predicted insert size before
                                   and after extension match.  Requires setting extend2.
                                   Enables the lowest possible false-positive rate.
        ecctadpole=f           (ecct) If reads fail to merge, error-correct with Tadpole
                                   and try again.  This happens prior to extend2.
        reassemble=t           If ecct is enabled, use Tadpole's reassemble mode for
                                   error correction.  Alternatives are pincer and tail.
        removedeadends         (shave) Remove kmers leading to dead ends.
        removebubbles          (rinse) Remove kmers in error bubbles.
        mindepthseed=3         (mds) Minimum kmer depth to begin extension.
        mindepthextend=2       (mde) Minimum kmer depth continue extension.
        branchmult1=20         Min ratio of 1st to 2nd-greatest path depth at high depth.
        branchmult2=3          Min ratio of 1st to 2nd-greatest path depth at low depth.
        branchlower=3          Max value of 2nd-greatest path depth to be considered low.
        ibb=t                  Ignore backward branches when extending.
        extra=<file>           A file or comma-delimited list of files of reads to use
                                   for kmer counting, but not for merging or output.
        prealloc=f             Pre-allocate memory rather than dynamically growing;
                                   faster and more memory-efficient for large datasets.
                                   A float fraction (0-1) may be specified, default 1.
        prefilter=0            If set to a positive integer, use a countmin sketch to
                                   ignore kmers with depth of that value or lower, to
                                   reduce memory usage.
        filtermem=0            Allows manually specifying prefilter memory in bytes, for
                                   deterministic runs.  0 will set it automatically.
        minprob=0.5            Ignore kmers with overall probability of correctness
                                   below this, to reduce memory usage.
        minapproxoverlap=26    For rem mode, do not merge reads if the extended reads
                                   indicate that the raw reads should have overlapped by
                                   at least this much, but no overlap was found.


    Bloom Filter Parameters (for kmer operations with less memory than Tadpole)
    *Note: These require more memory and should be run with bbmerge-auto.sh.*
        eccbloom=f            (eccb) If reads fail to merge, error-correct with bbcms
                                  and try again.
        testmerge=f           Test kmer counts around the read merge junctions.  If
                                  it appears that the merge created new errors, undo it.
                                  This reduces the false-positive rate, but not as much as
                                  rem or rsem.

    Java Parameters:
    -Xmx                     This will set Java's memory usage,
                                 overriding autodetection.
                                 For example, -Xmx400m will specify 400 MB RAM.
    -eoom                    This flag will cause the process to exit if an
                                 out-of-memory exception occurs.  Requires Java 8u92+.
    -da                      Disable assertions.

    Parameters
    ----------
    inp_fq_one_fnm : str
        Primary input
    inp_fq_two_fnm : str
        A second file
    outp_fnm : str
        File for merged reads (default: 'merged.fastq')
    outpu1_fnm : str
        File for unmerged reads (default: 'unmerged1.fastq')
    outpu2_fnm : str
        A second file (default: 'unmerged2.fastq')
    strictness : str
        False positive and merging rate (default: 'vloose=t')
    qin : int
        ASCII offset for input quality.  May be 33 (Sanger), 64
        (Illumina), or auto. (default: 33)
    args : tuple
        Values to be added as is to the command, each separated by a
        space from each other and the command
    kwargs : dict
        Key-value pairs to be added to the command, appending a double
        dash to the key, and replacing dash with underscore in the
        key, each separated by a space from each other and the
        command. Use a kwarg to provide trusted or untrusted contigs

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define BBMerge command
    command = " ".join(
        [
            "bbmerge.sh",
            f"in1={inp_fq_one_fnm}",
            f"in2={inp_fq_two_fnm}",
            f"out={outp_fnm}",
            f"outu1={outpu1_fnm}",
            f"outu2={outpu2_fnm}",
            strictness,
            f"qin={qin}",
        ]
    )
    for arg in args:
        command = " ".join([command, arg])
    for key, val in kwargs.items():
        command = " ".join(
            [
                command,
                "--" + key.replace("_", "-"),
                str(val),
            ]
        )

    # Run the command in the Docker image
    image = "ralatsdio/bbtools:v38.90"
    client = docker.from_env()
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )


def BBSplit(inp_fq_fNm, *args, **kwargs):
    """Using BBTools BBMap reformat.sh to split interleaved reads into
    two separated fastq files.

    BBSplit maps reads to multiple references simultaneously.  Outputs
    reads to a file for the reference they best match, with multiple
    options for dealing with ambiguous mappings.

    To index:     bbsplit.sh build=<1> ref_x=<reference fasta> ref_y=<another reference fasta>
    To map:       bbsplit.sh build=<1> in=<reads> out_x=<output file> out_y=<another output file>

    To be concise, and do everything in one command:
    bbsplit.sh ref=x.fa,y.fa in=reads.fq basename=o%.fq

    that is equivalent to
    bbsplit.sh build=1 in=reads.fq ref_x=x.fa ref_y=y.fa out_x=ox.fq out_y=oy.fq

    By default paired reads will yield interleaved output, but you can use the # symbol to produce twin output files.
    For example, basename=o%_#.fq will produce ox_1.fq, ox_2.fq, oy_1.fq, and oy_2.fq.


    Indexing Parameters (required when building the index):
    ref=<file,file>     A list of references, or directories containing fasta files.
    ref_<name>=<ref.fa> Alternate, longer way to specify references. e.g., ref_ecoli=ecoli.fa
                        These can also be comma-delimited lists of files; e.g., ref_a=a1.fa,a2.fa,a3.fa
    build=<1>           If multiple references are indexed in the same directory, each needs a unique build ID.
    path=<.>            Specify the location to write the index, if you don't want it in the current working directory.

    Input Parameters:
    build=<1>           Designate index to use.  Corresponds to the number specified when building the index.
    in=<reads.fq>       Primary reads input; required parameter.
    in2=<reads2.fq>     For paired reads in two files.
    qin=<auto>          Set to 33 or 64 to specify input quality value ASCII offset.
    interleaved=<auto>  True forces paired/interleaved input; false forces single-ended mapping.
                        If not specified, interleaved status will be autodetected from read names.

    Mapping Parameters:
    maxindel=<20>       Don't look for indels longer than this.  Lower is faster.  Set to >=100k for RNA-seq.
    minratio=<0.56>     Fraction of max alignment score required to keep a site.  Higher is faster.
    minhits=<1>         Minimum number of seed hits required for candidate sites.  Higher is faster.
    ambiguous=<best>    Set behavior on ambiguously-mapped reads (with multiple top-scoring mapping locations).
                           best   (use the first best site)
                           toss   (consider unmapped)
                           random   (select one top-scoring site randomly)
                           all   (retain all top-scoring sites.  Does not work yet with SAM output)
    ambiguous2=<best>   Set behavior only for reads that map ambiguously to multiple different references.
                        Normal 'ambiguous=' controls behavior on all ambiguous reads;
                        Ambiguous2 excludes reads that map ambiguously within a single reference.
                           best   (use the first best site)
                           toss   (consider unmapped)
                           all   (write a copy to the output for each reference to which it maps)
                           split   (write a copy to the AMBIGUOUS_ output for each reference to which it maps)
    qtrim=<true>        Quality-trim ends to Q5 before mapping.  Options are 'l' (left), 'r' (right), and 'lr' (both).
    untrim=<true>       Undo trimming after mapping.  Untrimmed bases will be soft-clipped in cigar strings.

    Output Parameters:
    out_<name>=<file>   Output reads that map to the reference <name> to <file>.
    basename=prefix%suffix     Equivalent to multiple out_%=prefix%suffix expressions, in which each % is replaced by the name of a reference file.
    bs=<file>           Write a shell script to 'file' that will turn the sam output into a sorted, indexed bam file.
    scafstats=<file>    Write statistics on how many reads mapped to which scaffold to this file.
    refstats=<file>     Write statistics on how many reads were assigned to which reference to this file.
                        Unmapped reads whose mate mapped to a reference are considered assigned and will be counted.
    nzo=t               Only print lines with nonzero coverage.

    ***** Notes *****
    Almost all BBMap parameters can be used; run bbmap.sh for more details.
    Exceptions include the 'nodisk' flag, which BBSplit does not support.
    BBSplit is recommended for fastq and fasta output, not for sam/bam output.
    When the reference sequences are shorter than read length, use Seal instead of BBSplit.

    Java Parameters:
    -Xmx                This will set Java's memory usage, overriding autodetection.
                        -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.
                        The max is typically 85% of physical memory.
    -eoom               This flag will cause the process to exit if an
                        out-of-memory exception occurs.  Requires Java 8u92+.
    -da                 Disable assertions.

    Parameters
    ----------
    inp_fq_fnm : str
        File containing interleaved reads in FASTQ format

    Returns
    -------
    str
        The container logs, either STDOUT, STDERR, or both, depending
        on the value of the stdout and stderr arguments

    """
    # Define BBMerge command
    basename = inp_fq_fnm.split(".", 1)
    outp1_fnm = f"{basename[0]}_R1.{basename[1]}"
    outp2_fnm = f"{basename[0]}_R2.{basename[1]}"
    command = " ".join(
        [
            "reformat.sh",
            f"in={inp_fq_fnm}",
            f"out1={outp1_fnm}",
            f"out2={outp2_fnm}",
        ]
    )
    for arg in args:
        command = " ".join([command, arg])
    for key, val in kwargs.items():
        command = " ".join(
            [
                command,
                "--" + key.replace("_", "-"),
                str(val),
            ]
        )

    # Run the command in the Docker image
    image = "ralatsdio/bbtools:v38.90"
    client = docker.from_env()
    hosting_dir = os.getcwd()
    working_dir = "/data"
    volumes = {hosting_dir: {"bind": working_dir, "mode": "rw"}}
    return client.containers.run(
        image,
        command=command,
        volumes=volumes,
        working_dir=working_dir,
    )
