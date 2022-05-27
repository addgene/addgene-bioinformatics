import contextlib
import glob
import gzip
import logging
import os
import shutil
import subprocess
import time

from Bio import Align, SeqIO
from Bio.Seq import Seq
import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import RunUtilities as ru


SPADES_OUTPUT_DIR = "spades"
SPADES_OUTPUT_FNM = "contigs.fasta"

APC_OUTPUT_DIR = "apc"
APC_BASE_FNM = "apc"
APC_OUTPUT_FNM = APC_BASE_FNM + ".1.fa"

# Logging configuration
root = logging.getLogger()
if not root.handlers:
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    root.addHandler(ch)

logger = logging.getLogger("AssemblyUtilities")
logger.setLevel(logging.INFO)


# TODO: Use pathlib?


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
    try:
        old_dir = os.getcwd()
        os.chdir(new_dir)
        yield
    finally:
        os.chdir(old_dir)


@contextlib.contextmanager
def timing(message):
    """Time a chunk of code and print timing messages.

    Parameters
    ----------
    message : str
        Message to log

    Returns
    -------
    None

    """
    try:
        start_time = time.time()
        logger.info(message + " started")
        yield
    finally:
        logger.info(message + f" finished in {time.time() - start_time} s")


def copy_reads(sequencing_data_dir, plate, well, working_dir, force=False):
    """Copy read files for a given plate and well to a working
    directory, if they don't already exist there.

    Parameters
    ----------
    sequencing_data_dir : str
        Path to directory containing sequencing data
    plate : str
        Plate identifier
    well : str
        Well identifier
    working_dir : str
        Path to working directory
    force : boolean
        Flag to force copying if read files exist

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
        if not os.path.exists(rd_fnm) or force:
            shutil.copy(rd_file_path, working_dir)

    return rd_fnms


def copy_adapters(sequencing_data_dir, working_dir, force):
    """Copy adapters file to a working directory, if it doesn't
    already exist there.

    Parameters
    ----------
    sequencing_data_dir : str
        Path to directory containing sequencing data
    working_dir : str
        Path to working directory
    force : boolean
        Flag to force copying if adapters file exists

    Returns
    -------
    adapter_fnm : str
        Adapter file name

    """
    # Copy each file, if needed
    adapter_fnm = "adapters.fa"
    adapter_file_path = os.path.join(sequencing_data_dir, adapter_fnm)
    if not os.path.exists(adapter_fnm) or force:
        shutil.copy(adapter_file_path, working_dir)

    return adapter_fnm


def preprocess_reads_using_bbtools(working_dir, rd1_fnm, rd2_fnm, force=False):
    """Preprocess reads uings BBTools BBduk, BBnorm, and BBMerge.

    Parameters
    ----------
    working_dir : str
        Directory in which preprocessing occurs
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    force : boolean
        Flag to force preprocessing if output files exist

    Returns
    -------
    rd1_unmerged_fnm : str
        File name of reads one, after optional preprocessing
    rd2_unmerged_fnm : str
        File name of reads two, after optional preprocessing
    rds_merged_fnm : str
        File name of Merged reads, after optional preprocessing

    """
    # Preprocess only if the output directory does not exist
    rd1_trimmed_fnm = rd1_fnm.replace(".fastq.gz", "_trimmed.fastq.gz")
    rd2_trimmed_fnm = rd2_fnm.replace(".fastq.gz", "_trimmed.fastq.gz")
    rd1_normalized_fnm = rd1_fnm.replace(".fastq.gz", "_normalized.fastq.gz")
    rd2_normalized_fnm = rd2_fnm.replace(".fastq.gz", "_normalized.fastq.gz")
    rd1_unmerged_fnm = rd1_fnm.replace(".fastq.gz", "_unmerged.fastq.gz")
    rd2_unmerged_fnm = rd2_fnm.replace(".fastq.gz", "_unmerged.fastq.gz")
    rds_merged_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
        ".fastq.gz", "_merged.fastq.gz"
    )
    with pushd(working_dir):

        # Run BBDuk, if needed
        if (
            not os.path.exists(rd1_trimmed_fnm)
            or not os.path.exists(rd2_trimmed_fnm)
            or force
        ):
            with timing("Trimming using BBDuk"):
                ru.BBDuk(
                    rd1_fnm,
                    rd2_fnm,
                    outp_fnm=rd1_trimmed_fnm,
                    outp2_fnm=rd2_trimmed_fnm,
                )

        # Run BBNorm, if needed
        if (
            not os.path.exists(rd1_normalized_fnm)
            or not os.path.exists(rd2_normalized_fnm)
            or force
        ):
            with timing("Normalizing using BBNorm"):
                ru.BBNorm(
                    rd1_trimmed_fnm,
                    rd2_trimmed_fnm,
                    outp_fnm=rd1_normalized_fnm,
                    outp2_fnm=rd2_normalized_fnm,
                )

        # Run BBMerge, if needed
        if (
            not os.path.exists(rds_merged_fnm)
            or not os.path.exists(rd1_unmerged_fnm)
            or not os.path.exists(rd2_unmerged_fnm)
            or force
        ):
            with timing("Merging using BBMerge"):
                ru.BBMerge(
                    rd1_normalized_fnm,
                    rd2_normalized_fnm,
                    outp_fnm=rds_merged_fnm,
                    outpu1_fnm=rd1_unmerged_fnm,
                    outpu2_fnm=rd2_unmerged_fnm,
                )

    return (
        rd1_trimmed_fnm,
        rd2_trimmed_fnm,
        rd1_normalized_fnm,
        rd2_normalized_fnm,
        rd1_unmerged_fnm,
        rd2_unmerged_fnm,
        rds_merged_fnm,
    )


def assemble_using_spades(working_dir, rd1_fnm, rd2_fnm, preprocess=True, force=False):
    """Assemble reads using SPAdes, then circularize using apc.

    Parameters
    ----------
    working_dir : str
        Directory in which the assembly occurs
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    force=False : boolean
        Flag to force assembly if output files exist

    Returns
    -------
    spades_seq :  Bio.Seq.Seq
        Sequence object returned by SPAdes
    apc_seq :  Bio.Seq.Seq
        Sequence object returned by apc
    rd1_unmerged_fnm : str
        File name of reads one, after optional preprocessing
    rd2_unmerged_fnm : str
        File name of reads two, after optional preprocessing
    rds_merged_fnm : str
        File name of Merged reads, after optional preprocessing

    """
    # Assemble only if the output directory does not exist
    spades_seq = None
    apc_seq = None
    rd1_unmerged_fnm = None
    rd2_unmerged_fnm = None
    rds_merged_fnm = None
    spades_dir = os.path.join(working_dir, SPADES_OUTPUT_DIR)
    output_pth = os.path.join(SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
    with pushd(working_dir):

        # Optionally preprocess, if needed
        if preprocess:
            with timing("Preprocessing using BBTools"):
                (
                    rd1_trimmed_fnm,
                    rd2_trimmed_fnm,
                    rd1_normalized_fnm,
                    rd2_normalized_fnm,
                    rd1_unmerged_fnm,
                    rd2_unmerged_fnm,
                    rds_merged_fnm,
                ) = preprocess_reads_using_bbtools(
                    working_dir, rd1_fnm, rd2_fnm, force=force
                )

            # Run SPAdes with merged reads, if needed
            if not os.path.exists(output_pth) or force:
                with timing("Assembling preprocessed reads using SPAdes"):
                    ru.spades(
                        rd1_unmerged_fnm,
                        rd2_unmerged_fnm,
                        SPADES_OUTPUT_DIR,
                        merged=rds_merged_fnm,
                        cov_cutoff=100,
                    )

        else:
            # Run SPAdes without merged reads, if needed
            if not os.path.exists(output_pth) or force:
                with timing("Assembling original reads using SPAdes"):
                    ru.spades(
                        rd1_fnm,
                        rd2_fnm,
                        SPADES_OUTPUT_DIR,
                        cov_cutoff=100,
                    )

        # Copy SPAdes output file, and parse
        if os.path.exists(output_pth):
            shutil.copy(output_pth, ".")
            spades_seq = [
                seq_rcd for seq_rcd in SeqIO.parse(SPADES_OUTPUT_FNM, "fasta")
            ][0].seq

        # Run apc
        with timing("Circularize using apc"):
            ru.apc(APC_BASE_FNM, SPADES_OUTPUT_FNM)
        logger.info(
            f"Results of assembling usng SPAdes in directory: {os.getcwd()}"
        )

        # Parse apc output file
        if os.path.exists(APC_OUTPUT_FNM):
            apc_seq = [seq_rcd for seq_rcd in SeqIO.parse(APC_OUTPUT_FNM, "fasta")][
                0
            ].seq

        return (
            spades_seq,
            apc_seq,
            rd1_unmerged_fnm,
            rd2_unmerged_fnm,
            rds_merged_fnm,
        )


def compute_coverage_using_bbtools(working_dir, rd1_fnm, rd2_fnm, ref_fnm, force=False):
    """Compute coverage uings BBTools BBMap.

    Parameters
    ----------
    working_dir : str
        Directory in which preprocessing occurs
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    ref_fnm : str
        Name of reference sequence file
    force=False : boolean
        Flag to force computing if output files exist

    Returns
    -------
    coverage : float
        Average number of reads that map to the reference sequence
    out_fnm : str
        File name of reads mapped to the reference sequence
    covstats_fnm : str
        File name of the file to which coverage statistics are written

    """
    # Compute coverage, if needed
    coverage = 0
    out_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
        ".fastq.gz", "_output.fastq.gz"
    )
    covstats_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
        ".fastq.gz", "_covstats.txt"
    )
    with pushd(working_dir):
        if not os.path.exists(out_fnm) or not os.path.exists(covstats_fnm) or force:

            # Run BBMap
            with timing("Mapping using BBMap"):
                ru.BBMap(rd1_fnm, rd2_fnm, out_fnm, ref_fnm, covstats_fnm)

            # Parse BBMap coverage statistics
            cp = subprocess.run(
                f"head -n 2 {covstats_fnm} | tail -n 1 | cut -wf 2",
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            coverage = float(cp.stdout.strip())

    return coverage, out_fnm, covstats_fnm


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
    k_mers : dct
        Emtpy, or returned by this method

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
        k_mer_rc = Seq(k_mer).reverse_complement()
        if k_mer not in k_mers and k_mer_rc not in k_mers:
            k_mers[k_mer] = {}
            k_mers[k_mer]["src"] = set([seq_id])
            k_mers[k_mer]["cnt"] = 1
        elif k_mer in k_mers:
            k_mers[k_mer]["src"].add(seq_id)
            k_mers[k_mer]["cnt"] += 1
        elif k_mer_rc in k_mers:
            k_mers[k_mer_rc]["src"].add(seq_id)
            k_mers[k_mer_rc]["cnt"] += 1
        else:
            raise Exception("Should not get here")

    return k_mers


def count_k_mers_in_rds(rd_fnm, k_mer_len=25, k_mers=None, seq_rcds=None):
    """Count k-mers of a specified length in a gzipped file of reads
    in the specified format. Optionally update dictionary and list
    returned by this method.

    Parameters
    ----------
    rd_fnm : str
        The name of the gzipped read file
    k_mers : dict
        Empty, or returned by count_k_mers_in_seq()
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Empty, or returned by count_k_mers_in_rds()
    k_mer_len : int
        The length of k-mers to count

    Returns
    -------
    k_mers : dict
        Containing k-mer keys and count and source sequence record
        index values
    seq_rcds : list(Bio.SeqRecord.SeqRecord)
        Contains sequence records in which k-mers were counted

    """
    if k_mers is None:
        k_mers = {}
    if seq_rcds is None:
        seq_rcds = []
    read_format, is_gzipped = ru.get_bio_read_format(rd_fnm)
    _open = open
    if is_gzipped:
        _open = gzip.open
    with _open(rd_fnm, "rt") as f:
        seq_rcd_gen = SeqIO.parse(f, format=read_format)
        i_seq = -1
        for seq_rcd in seq_rcd_gen:
            i_seq += 1
            k_mers = count_k_mers_in_seq(seq_rcd.seq, seq_id=i_seq, k_mers=k_mers)
            seq_rcds.append(seq_rcd)

    return k_mers, seq_rcds


def write_k_mer_counts_in_rds(k_mers_in_rd1, k_mers_in_rd2, k_mer_counts_fnm):
    """Write the k-mer counts corresponding to paired reads to the
    specified file.

    Parameters
    ----------
    k_mers_in_rd1 : dct
        Dictionary returned by count_k_mers_in_seq()
    k_mers_in_rd2 : dct
        Dictionary returned by count_k_mers_in_seq()
    k_mer_counts_fnm
        File name to which to write k-mers and their counts

    Returns
    -------
    None

    """
    with open(k_mer_counts_fnm, "w") as f:
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


def read_k_mer_counts(k_mer_counts_fnm, seq_id=0, k_mers=None):
    """Read the k-mer counts corresponding to paired reads from the
    specified file. Optionally update dictionary returned by
    count_k_mers_in_seq().

    Parameters
    ----------
    k_mers_counts_fnm : str
        The name of the file containing k-mers and their counts
    seq_id : int
        The sequence identifer
    k_mers : None or dct
        None or dictionary returned by count_k_mers_in_seq()

    Returns
    -------
    k_mers : dict
        Dictionay containing k-mer keys and counts, and source
        sequence identifiers values
    """
    if k_mers is None:
        k_mers = {}
    with open(k_mer_counts_fnm) as f:
        for ln in f:
            flds = ln.split()
            k_mer = flds[0]
            cnt_rd1 = int(flds[1])
            cnt_rd2 = int(flds[2])
            if k_mer in k_mers:
                logger.ingof(f"Second occurance of k-mer: {k_mer} unexpected")
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
    rd1_wr_file_name : str
        File name containing reads one with repeats
    rd1_wo_file_name : str
        File name containing reads one without repeats
    rd2_wr_file_name : str
        File name containing reads two with repeats
    rd2_wo_file_name : str
        File name containing reads two without repeats

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
    return rd1_wr_file_name, rd1_wo_file_name, rd2_wr_file_name, rd2_wo_file_name


def find_seq_rcds_for_cnt(k_mer_cnt_rds, min_n_clusters=2, max_n_clusters=8):
    """TODO: Complete

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
    i_rcds : set(int)
        Index of sequence records in a cluster
    kmeans : sklearn.cluster.KMeans
        Object containing result of k-means estimation

    """
    # Cluster counts in the expected number of clusters
    opt_kmeans = None
    opt_s_score = -1
    opt_n_clusters = min_n_clusters
    # X = (k_mer_cnt_rds / coverage).reshape(-1, 1)
    X = k_mer_cnt_rds.reshape(-1, 1)
    with timing("Selecting number of clusters"):
        for n_clusters in range(min_n_clusters, max_n_clusters + 1):
            with timing(f"Fitting and scoring {n_clusters} clusters"):
                kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(X)
                s_score = silhouette_score(X, kmeans.labels_, random_state=0)
                logger.info(
                    f"With {n_clusters} clusters found silhouette score {s_score}"
                )
                if s_score > opt_s_score:
                    opt_kmeans = kmeans
                    opt_s_score = s_score
                    opt_n_clusters = n_clusters

    logger.info(
        f"Optimum number of clusters {opt_n_clusters} has silhouette score {opt_s_score}"
    )

    # Identify the cluster and sequence records (reads) corresponding
    # to the expected count
    """
    label = kmeans.predict(np.array(K_MER_CNT_REP).reshape(-1, 1))[0]
    i_rcds = set()
    keys = np.array(list(k_mers.keys()))
    for key in keys[kmeans.labels_ == label]:
        i_rcds = i_rcds.union(k_mers[key]["src"])

    return i_rcds, kmeans
    """
    return opt_kmeans, opt_s_score, opt_n_clusters


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
    read_wr_file_name : str
        File name of reads with repeats
    read_wo_file_name : str
        File name of reads without repeats

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
def ghi(initial_seq, rd1_file_name, rd2_file_name):
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]
    OUTER_DISTANCE = 500
    FRAGMENT_LEN = OUTER_DISTANCE

    # Count k-mers in the random sequence, doubled to represent
    # circular DNA
    # TODO: Fix
    with timing("Counting k_mers in initial sequence"):
        k_mers_seq = count_k_mers_in_seq(initial_seq + initial_seq)

    # Count k-mers in the paired reads, and write the result to a
    # file
    with timing("Counting k_mers in reads"):
        k_mers_rd1, seq_rcds_rd1 = count_k_mers_in_rds(rd1_file_name)
        k_mers_rd2, seq_rcds_rd2 = count_k_mers_in_rds(rd2_file_name)
    with timing("Writing k_mers and counts in reads"):
        write_k_mer_counts_in_rds(k_mers_rd1, k_mers_rd2, BASE_FILE_NAME + "_cnt.txt")

    # Collect k-mer counts, and compute coverage
    k_mer_cnt_seq = collect_k_mer_cnt(k_mers_seq, seq_cnt=2)
    k_mer_cnt_rd1 = collect_k_mer_cnt(k_mers_rd1, seq_cnt=1)
    k_mer_cnt_rd2 = collect_k_mer_cnt(k_mers_rd2, seq_cnt=1)
    coverage_rd1 = int(np.sum(k_mer_cnt_rd1) / np.sum(k_mer_cnt_seq))
    coverage_rd2 = int(np.sum(k_mer_cnt_rd2) / np.sum(k_mer_cnt_seq))

    # Separate reads into those containing a k-mer with the
    # expected count, and all others
    with timing("Writing reads based on read count"):
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

    # Assemble paired reads with repeats using SSAKE
    with timing("Assembling paired reads with repeats using SSAKE"):
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

    # Assemble paired reads with trusted contigs using SPAdes
    with timing(
        "Assembling paired reads repeats and with trusted contigs using SPAdes"
    ):
        spades_wc_out_dir = BASE_FILE_NAME + "_spades_wc"
        trusted_contigs_fnm = os.path.join(
            ssake_out_base_name, ssake_out_base_name + "_scaffolds.fa"
        )
        if os.path.exists(spades_wc_out_dir):
            shutil.rmtree(spades_wc_out_dir)
        ru.spades(
            rd1_file_name,
            rd2_file_name,
            spades_wc_out_dir,
            trusted_contigs_fnm=trusted_contigs_fnm,
        )


def align_assembly_output(aligner, working_dir, seq):
    """Align SPAdes and apc assembly output to a sequence.

    Parameters
    ----------
    aligner : Align.PairwiseAligner
        A Biopython pairwise aligner
    working_dir
        Directory containing output of an assembly case
    seq : Bio.Seq.Seq
        A Biopython sequence to which to align

    Returns
    -------
    spd_scr
        Alignment score for the SPAdes output
    apc_scr
        Alignment score for the apc output

    """
    case = os.path.basename(working_dir)
    with timing(f"Aligning {case} assemblies"):

        # Read the multi-FASTA SPAdes output file, and align
        output_pth = os.path.join(working_dir, SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
        spd_scr = -1.0
        if os.path.exists(output_pth):
            spd_scr = 0.0
            seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
            if len(seq_rcds) > 0:
                spd_scr = aligner.score(seq + seq, seq_rcds[0].seq)

        # Read the FASTA apc output file, and align
        output_pth = os.path.join(working_dir, APC_OUTPUT_FNM)
        apc_scr = -1.0
        if os.path.exists(output_pth):
            apc_scr = 0.0
            seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
            if len(seq_rcds) > 0:
                apc_scr = aligner.score(seq + seq, seq_rcds[0].seq)

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
