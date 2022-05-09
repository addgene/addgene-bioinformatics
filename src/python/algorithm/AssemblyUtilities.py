import contextlib
import glob
import gzip
import os
import shutil
import time

from Bio import Align, SeqIO
import numpy as np
from sklearn.cluster import KMeans

import RunUtilities as ru


SPADES_OUTPUT_DIR = "spades"
SPADES_OUTPUT_FNM = "contigs.fasta"

APC_OUTPUT_DIR = "apc"
APC_BASE_FNM = "apc"
APC_OUTPUT_FNM = APC_BASE_FNM + ".1.fa"


@contextlib.contextmanager
def pushd(new_dir):
    """Change to a new working directory, and back.

    Parameters
    ----------
    new_dir : str
        Path of the target directory

    Returns
    -------
    None

    """
    old_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(old_dir)


def copy_actual_reads(plate, well, working_dir, sequencing_data_dir):
    """Copy actual read files for a given plate and well to a working
    directory, if they don't already exist there.

    Parameters
    ----------
    plate : str
        Plate identifier
    well : str
        Well identifier
    working_dir : str
        Path to working directory
    sequencing_data_dir : str
        Path to directory containing sequencing data

    Returns
    -------
    rd_fnms : tuple(str)
        Tuple containing read file names

    """
    # Find the plate directory given that we only know the Addgene,
    # but not the seqWell, identifier
    plate_dir = os.path.basename(
        glob.glob(os.path.join(sequencing_data_dir, plate + "*"))[0]
    )

    # Follow the read file naming convention
    rd1_fnm = plate_dir.replace("_FASTQ", "") + "_" + well + "_R1_001.fastq.gz"
    rd2_fnm = rd1_fnm.replace("R1", "R2")

    # Copy each file, if needed
    rd_fnms = (rd1_fnm, rd2_fnm)
    for rd_fnm in rd_fnms:
        rd_file_path = os.path.join(sequencing_data_dir, plate_dir, rd_fnm)
        if not os.path.exists(rd_fnm):
            shutil.copy(rd_file_path, working_dir)

    return rd_fnms


def assemble_using_spades(working_diur, rd1_fnm, rd2_fnm, force=False):
    """Assemble reads using SPAdes, then circularize using apc.

    Parameters
    ----------
    working_diur : str
        Directory in which the assembly occurs
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    force=False : boolean
        Flag to force assembly if output directory exists

    Returns
    -------
    None

    """
    # Assemble only if the output directory does not exist
    spades_dir = os.path.join(working_diur, SPADES_OUTPUT_DIR)
    if not os.path.exists(spades_dir) or force:
        with pushd(working_diur):

            # Run SPAdes
            # TODO: Check parameters
            start_time = time.time()
            print("Assembling using SPAdes ...", end=" ", flush=True)
            command = ru.spades(rd1_fnm, rd2_fnm, SPADES_OUTPUT_DIR, cov_cutoff=100)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)

            # Copy SPAdes output file
            output_pth = os.path.join(SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
            if os.path.exists(output_pth):
                shutil.copy(output_pth, ".")

            # Run apc
            start_time = time.time()
            print("Circularize using apc ...", end=" ", flush=True)
            command = ru.apc(APC_BASE_FNM, SPADES_OUTPUT_FNM)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)
            print("Directory: {0}".format(os.getcwd()))


def count_k_mers_in_seq(seq, seq_id=0, k_mer_len=25, k_mers=None):
    """Count k-mers of a specified length in a sequence with a
    specified identifier. Optionally update dictionary returned by
    this method.

    Parameters
    ----------
    seq : Bio.Seq.Seq
        The sequence in which to count
    seq_id : int
        The source sequence identifer
    k_mer_len : int
        The length of k-mers to count
    k_mers : None or dct
        None or dictionary returned by this method

    Returns
    -------
    k_mers : dct
        Dictionay containing k-mer keys and counts and source sequence
        identifier values

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
    in the specified format. Optionally update dictionary and list
    returned by this method.

    Parameters
    ----------
    rd_fNm : str
        The name of the gzipped read file
    k_mer_len : int
        The length of k-mers to count
    k_mers : None or dict
        None or dictionary returned by count_k_mers_in_seq()
    seq_rcds : None or list(Bio.SeqRecord.SeqRecord)
        None or list of sequence records in which k-mers were counted

    Returns
    -------
    tuple(dict, list(Bio.SeqRecord.SeqRecord))
        Dictionay containing k-mer keys and count and source sequence
        record index values, and list of sequence records in which
        k-mers were counted
    """
    if k_mers is None:
        k_mers = {}
    if seq_rcds is None:
        seq_rcds = []
    read_format, is_gzipped = ru.get_bio_read_format(rd_fNm)
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
    """Write the k-mer counts corresponding to paired reads to the
    specified file.

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
                "{0} {1:10d} {2:10d}\n".format(
                    k_mer, int(k_mer_cnt_rd1), int(k_mer_cnt_rd2)
                )
            )


def read_k_mer_counts(k_mer_counts_fNm, seq_id=0, k_mers=None):
    """Read the k-mer counts corresponding to paired reads from the
    specified file. Optionally update dictionary returned by
    count_k_mers_in_seq().

    Parameters
    ----------
    k_mers_counts_fNm : str
        The name of the file containing k-mers and their counts
    seq_id : int
        The sequence identifer
    k_mers : None or dct
        None or dictionary returned by count_k_mers_in_seq()

    Returns
    -------
    dict
        Dictionay containing k-mer keys and counts, and source
        sequence identifiers values
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


def collect_k_mer_cnt(k_mers, seq_cnt=2):
    """Collect k-mer counts, assuming that we counted in a doubled
    sequence. Set sequence count otherwise.

    Parameters
    ----------
    k_mers : dict
        Dictionary returned by count_k_mers_in_seq()
    seq_cnt : int
        Number of times the sequence was repeated for counting
        (default: 2)

    Returns
    -------
    numpy.ndarray
        Counts of k_mers

    """
    return np.array([val["cnt"] / seq_cnt for val in k_mers.values()])


def write_paired_reads_for_cnt(
    k_mer_cnt_rd1,
    coverage_rd1,
    k_mers_rd1,
    seq_rcds_rd1,
    k_mer_cnt_rd2,
    coverage_rd2,
    k_mers_rd2,
    seq_rcds_rd2,
):
    # TODO: Update.
    """Write paired reads for counts.

    Parameters
    ----------
    k_mer_cnt_rd1 : numpy.ndarray
        Counts of k_mers
    coverage_rd1 : int
        Read one coverage
    k_mers_rd1 : dict
        Contains k-mer keys and count and source sequence record index
        values
    seq_rcds_rd1 : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted
    k_mer_cnt_rd2  : numpy.ndarray
        Counts of k_mers
    coverage_rd2 : int
        Read two coverage
    k_mers_rd2 : dict
        Contains k-mer keys and count and source sequence record index
        values
    seq_rcds_rd2 : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted

    Returns
    -------
    tuple(str)
        Read one and two file names with and without repeats

    """
    # Find common sequence records
    if len(seq_rcds_rd1) != len(seq_rcds_rd2):
        raise (Exception("Number of sequence records for each pair must agree"))
    i_rcds_rd1, _ = find_seq_rcds_for_cnt(k_mer_cnt_rd1, coverage_rd1, k_mers_rd1)
    i_rcds_rd2, _ = find_seq_rcds_for_cnt(k_mer_cnt_rd2, coverage_rd2, k_mers_rd2)
    i_rcds = i_rcds_rd1.intersection(i_rcds_rd2)

    # Write read one and two files with and without repeats
    rd1_wr_file_name, rd1_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd1, "rd1"
    )
    rd2_wr_file_name, rd2_wo_file_name = write_reads_for_cnt(
        i_rcds, seq_rcds_rd2, "rd2"
    )
    return (rd1_wr_file_name, rd1_wo_file_name, rd2_wr_file_name, rd2_wo_file_name)


# TODO: Step through to validate
def find_seq_rcds_for_cnt(k_mer_cnt_rds, coverage, k_mers):
    """Use k-means clustering to identify reads containing a k-mer
    with the expected count.

    Parameters
    ----------
    k_mer_cnt_rds : numpy.ndarray
        Counts of k_mers
    coverage : int
        Read coverage
    k_mers : dict
        Contains k-mer keys and count and source sequence record index
        values

    Returns
    -------
    tuple(set(int), sklearn.cluster.KMeans)
        Index of sequence records in a cluster

    """
    # TODO: Settle
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]

    # Cluster counts in the expected number of clusters
    kmeans = KMeans(n_clusters=len(K_MER_CNT_REP) + 1, random_state=0).fit(
        (k_mer_cnt_rds / coverage).reshape(-1, 1)

    # Identify the cluster and sequence records (reads) corresponding
    # to the expected count
    label = kmeans.predict(np.array(K_MER_CNT_REP).reshape(-1, 1))[0]
    i_rcds = set()
    keys = np.array(list(k_mers.keys()))
    for key in keys[kmeans.labels_ == label]:
        i_rcds = i_rcds.union(k_mers[key]["src"])

    return i_rcds, kmeans


def write_reads_for_cnt(i_rcds, seq_rcds, case):
    """Write a FASTQ file containing the sequence records (reads)
    specified, and a file containing all other sequence records
    (reads).

    Parameters
    ----------
    i_rcds : set(int)
        Index of sequence records in a cluster
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records
    case : str
        Label for read files

    Returns
    -------
    tuple(str)
        Names of read files

    """
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    NUMBER_PAIRS = 25000

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


# TODO: Name and review
def ghi(rd1_file_name, rd2_file_name):
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]
    OUTER_DISTANCE = 500
    FRAGMENT_LEN = OUTER_DISTANCE

    # Count k-mers in the paired reads, and write the result to a
    # file
    start_time = time.time()
    print("Counting k_mers in reads ...", end=" ", flush=True)
    k_mers_rd1, seq_rcds_rd1 = count_k_mers_in_rds(rd1_file_name)
    k_mers_rd2, seq_rcds_rd2 = count_k_mers_in_rds(rd2_file_name)
    print("done in {0} s".format(time.time() - start_time), flush=True)
    start_time = time.time()
    print("Writing k_mers and counts in reads ...", end=" ", flush=True)
    write_k_mer_counts_in_rds(k_mers_rd1, k_mers_rd2, BASE_FILE_NAME + "_cnt.txt")
    print("done in {0} s".format(time.time() - start_time), flush=True)

    # Collect k-mer counts, and compute coverage
    # TODO: How to count kmers in the sequence before it is known
    # k_mer_cnt_seq = collect_k_mer_cnt(k_mers_seq, seq_cnt=2)
    k_mer_cnt_rd1 = collect_k_mer_cnt(k_mers_rd1, seq_cnt=1)
    k_mer_cnt_rd2 = collect_k_mer_cnt(k_mers_rd2, seq_cnt=1)
    # coverage_rd1 = int(np.sum(k_mer_cnt_rd1) / np.sum(k_mer_cnt_seq))
    # coverage_rd2 = int(np.sum(k_mer_cnt_rd2) / np.sum(k_mer_cnt_seq))
    # coverage_rd1 = 1
    # coverage_rd2 = 1

    # Separate reads into those containing a k-mer with the
    # expected count, and all others
    print("Writing reads based on read count ...", end=" ", flush=True)
    # TODO: What to do with this r_k_mers
    # r_k_mer_exp_rd1 = str(r_k_mers[K_MER_CNT_REP.index(EXP_CNT)])
    # r_k_mer_exp_rd2 = r_k_mer_exp_rd1[-1::-1]
    (
        rd1_wr_file_name,
        rd1_wo_file_name,
        rd2_wr_file_name,
        rd2_wo_file_name,
    ) = write_paired_reads_for_cnt(
        k_mer_cnt_rd1,
        coverage_rd1,
        k_mers_rd1,
        seq_rcds_rd1,
        k_mer_cnt_rd2,
        coverage_rd2,
        k_mers_rd2,
        seq_rcds_rd2,
        K_MER_CNT_REP,
        EXP_CNT,
    )
    print("done in {0} s".format(time.time() - start_time), flush=True)
    # print("r_k_mer_exp_rd1: {0}".format(r_k_mer_exp_rd1))
    # print("r_k_mer_exp_rd2: {0}".format(r_k_mer_exp_rd2))

    # Assemble paired reads with repeats using SSAKE
    start_time = time.time()
    print("Assembling paired reads with repeats using SSAKE ...", end=" ", flush=True)
    ssake_out_base_name = BASE_FILE_NAME + "_ssake"
    if os.path.exists(ssake_out_base_name):
        shutil.rmtree(ssake_out_base_name)
        ru.ssake(
            rd1_wr_file_name,
            rd2_wr_file_name,
            ssake_out_base_name,
            FRAGMENT_LEN,
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
        )  # -v
    print("done in {0} s".format(time.time() - start_time), flush=True)

    # Assemble paired reads with trusted contigs using SPAdes
    start_time = time.time()
    print(
        "Assembling paired reads repeats "
        + "and with trusted contigs using SPAdes ...",
        end=" ",
        flush=True,
    )
    spades_wc_out_dir = BASE_FILE_NAME + "_spades_wc"
    trusted_contigs_fNm = os.path.join(
        ssake_out_base_name, ssake_out_base_name + "_scaffolds.fa"
    )
    if os.path.exists(spades_wc_out_dir):
        shutil.rmtree(spades_wc_out_dir)
        command = ru.spades(
            rd1_file_name,
            rd2_file_name,
            spades_wc_out_dir,
            trusted_contigs_fNm=trusted_contigs_fNm,
        )
    print("done in {0} s".format(time.time() - start_time), flush=True)
    print("Command: {0}".format(command), flush=True)


def align_assembly_output(aligner, working_diur, seq):
    """Align SPAdes and apc assembly output to a sequence.

    Parameters
    ----------
    aligner : Align.PairwiseAligner
        A Biopython pairwise aligner
    working_diur
        Directory containing output of an assembly case
    seq : Bio.Seq.Seq
        A Biopython sequence to which to align

    Returns
    -------
    tuple(float)
        Alignment score for the SPAdes and apc output

    """
    # Read the multi-FASTA SPAdes output file, and align
    start_time = time.time()
    case = os.path.basename(working_diur)
    print("Aligning {0} assemblies ...".format(case), end=" ", flush=True)
    output_pth = os.path.join(working_diur, SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
    spd_scr = -1.0
    if os.path.exists(output_pth):
        spd_scr = 0.0
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
        if len(seq_rcds) > 0:
            spd_scr = aligner.score(seq + seq, seq_rcds[0].seq)

    # Read the FASTA apc output file, and align
    output_pth = os.path.join(working_diur, APC_OUTPUT_FNM)
    apc_scr = -1.0
    if os.path.exists(output_pth):
        apc_scr = 0.0
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
        if len(seq_rcds) > 0:
            apc_scr = aligner.score(seq + seq, seq_rcds[0].seq)
    print("done in {0} s".format(time.time() - start_time), flush=True)
    return spd_scr, apc_scr


def create_aligner(config):
    """Creates a pairwise aligner with the specified configuration.

    Parameters
    ----------
    config : dict
        A pairwise aligner configuration

    Returns
    -------
    aligner : Align.PairwiseAligner
        A Biopython pairwise aligner

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
