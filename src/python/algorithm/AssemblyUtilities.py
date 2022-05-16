import contextlib
import glob
import gzip
import os
import shutil
import subprocess
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
    old_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(old_dir)


# TODO: Add context manager for timing?


def copy_reads(sequencing_data_dir, plate, well, working_dir):
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


def copy_adapters(sequencing_data_dir, working_dir):
    """Copy adapters file to a working directory, if it doesn't
    already exist there.

    Parameters
    ----------
    sequencing_data_dir : str
        Path to directory containing sequencing data
    working_dir : str
        Path to working directory

    Returns
    -------
    adapter_fnm : str
        Adapte file name

    """
    # Copy each file, if needed
    adapter_fnm = "adapters.fa"
    adapter_file_path = os.path.join(sequencing_data_dir, adapter_fnm)
    if not os.path.exists(adapter_fnm):
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
    force=False : boolean
        Flag to force preprocessing if output directory exists

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
    if not os.path.exists(working_dir) or force:
        with pushd(working_dir):

            # Run BBDuk
            start_time = time.time()
            print("Trimming using BBDuk ...", end=" ", flush=True)
            rd1_trimmed_fnm = rd1_fnm.replace(".fastq.gz", "_trimmed.fastq.gz")
            rd2_trimmed_fnm = rd2_fnm.replace(".fastq.gz", "_trimmed.fastq.gz")
            ru.BBDuk(
                rd1_fnm, rd2_fnm, outp_fnm=rd1_trimmed_fnm, outp2_fnm=rd2_trimmed_fnm
            )
            print("done in {0} s".format(time.time() - start_time), flush=True)

            # Run BBNorm
            start_time = time.time()
            print("Normalizing using BBNorm ...", end=" ", flush=True)
            rd1_normalized_fnm = rd1_fnm.replace(".fastq.gz", "_normalized.fastq.gz")
            rd2_normalized_fnm = rd2_fnm.replace(".fastq.gz", "_normalized.fastq.gz")
            ru.BBNorm(
                rd1_trimmed_fnm,
                rd2_trimmed_fnm,
                outp_fnm=rd1_normalized_fnm,
                outp2_fnm=rd2_normalized_fnm,
            )
            print("done in {0} s".format(time.time() - start_time), flush=True)

            # Run BBMerge
            start_time = time.time()
            print("Merging using BBMerge ...", end=" ", flush=True)
            rds_merged_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
                ".fastq.gz", "_merged.fastq.gz"
            )
            rd1_unmerged_fnm = rd1_fnm.replace(".fastq.gz", "_unmerged.fastq.gz")
            rd2_unmerged_fnm = rd2_fnm.replace(".fastq.gz", "_unmerged.fastq.gz")
            ru.BBMerge(
                rd1_normalized_fnm,
                rd2_normalized_fnm,
                outp_fnm=rds_merged_fnm,
                outpu1_fnm=rd1_unmerged_fnm,
                outpu2_fnm=rd2_unmerged_fnm,
            )
            print("done in {0} s".format(time.time() - start_time), flush=True)

        return rd1_unmerged_fnm, rd2_unmerged_fnm, rds_merged_fnm

    else:

        return ()


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
        Flag to force assembly if output directory exists

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
    spades_dir = os.path.join(working_dir, SPADES_OUTPUT_DIR)
    if not os.path.exists(spades_dir) or force:
        with pushd(working_dir):

            # Optionally preprocess
            rds_merged_fnm = None
            rd1_unmerged_fnm = None
            rd2_unmerged_fnm = None
            if preprocess:
                (
                    rds_merged_fnm,
                    rd1_unmerged_fnm,
                    rd2_unmerged_fnm,
                ) = preprocess_reads_using_bbtools(
                    working_dir, rd1_fnm, rd2_fnm, force=force
                )

                # Run SPAdes with merged reads
                start_time = time.time()
                print("Assembling using SPAdes ...", end=" ", flush=True)
                command = ru.spades(
                    rd1_unmerged_fnm,
                    rd2_unmerged_fnm,
                    SPADES_OUTPUT_DIR,
                    merged=rds_merged_fnm,
                    cov_cutoff=100,
                )
                print("done in {0} s".format(time.time() - start_time), flush=True)

            else:
                # Run SPAdes without merged reads
                start_time = time.time()
                print("Assembling using SPAdes ...", end=" ", flush=True)
                command = ru.spades(
                    rd1_fnm,
                    rd2_fnm,
                    SPADES_OUTPUT_DIR,
                    cov_cutoff=100,
                )
                print("done in {0} s".format(time.time() - start_time), flush=True)

            # Copy SPAdes output file, and parse
            spades_seq = None
            output_pth = os.path.join(SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
            if os.path.exists(output_pth):
                shutil.copy(output_pth, ".")
                spades_seq = [
                    seq_rcd for seq_rcd in SeqIO.parse(SPADES_OUTPUT_FNM, "fasta")
                ][0].seq

            # Run apc
            start_time = time.time()
            print("Circularize using apc ...", end=" ", flush=True)
            command = ru.apc(APC_BASE_FNM, SPADES_OUTPUT_FNM)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Directory: {0}".format(os.getcwd()))

            # Parse apc output file
            apc_seq = None
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
        Flag to force preprocessing if output directory exists

    Returns
    -------
    coverage : float
    out_fnm : str
        File name of reads mapped to the reference sequence
    covstats_fnm : str
        File name of the file to which coverage statistics are written

    """
    # Preprocess only if the output directory does not exist
    if not os.path.exists(working_dir) or force:
        with pushd(working_dir):

            # Run BBMap
            start_time = time.time()
            print("Mapping using BBMap ...", end=" ", flush=True)
            out_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
                ".fastq.gz", "_output.fastq.gz"
            )
            covstats_fnm = rd1_fnm.replace("_R1", "_Rs").replace(
                ".fastq.gz", "_covstats.txt"
            )
            ru.BBMap(rd1_fnm, rd2_fnm, out_fnm, ref_fnm, covstats_fnm)
            print("done in {0} s".format(time.time() - start_time), flush=True)

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

    else:

        return ()


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


def count_k_mers_in_rds(rd_fnm, k_mer_len=25, k_mers=None, seq_rcds=None):
    """Count k-mers of a specified length in a gzipped file of reads
    in the specified format. Optionally update dictionary and list
    returned by this method.

    Parameters
    ----------
    rd_fnm : str
        The name of the gzipped read file
    k_mer_len : int
        The length of k-mers to count
    k_mers : None or dict
        None or dictionary returned by count_k_mers_in_seq()
    seq_rcds : None or list(Bio.SeqRecord.SeqRecord)
        None or list of sequence records in which k-mers were counted

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
            k_mers = count_k_mers_in_seq(seq_rcd.seq, i_seq, k_mers=k_mers)
            seq_rcds.append(seq_rcd)

    return k_mers, seq_rcds


def write_k_mer_counts_in_rds(k_mers_in_rd1, k_mers_in_rd2, k_mer_counts_fnm):
    """Write the k-mer counts corresponding to paired reads to the
    specified file.

    Parameters
    ----------
    k_mers_in_rd1 : dct
        dictionary returned by count_k_mers_in_seq()
    k_mers_in_rd2 : dct
        dictionary returned by count_k_mers_in_seq()
    k_mer_counts_fnm
        file name to which to write k-mers and their counts

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


# TODO: Step through to validate
def find_seq_rcds_for_cnt(k_mer_cnt_rds, coverage, k_mers, n_clusters):
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
    i_rcds : set(int)
        Index of sequence records in a cluster
    kmeans : sklearn.cluster.KMeans
        Object containing result of k-means estimation

    """
    # TODO: Settle
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]

    # Cluster counts in the expected number of clusters
    kmeans = KMeans(n_clusters=len(K_MER_CNT_REP) + 1, random_state=0).fit(
        (k_mer_cnt_rds / coverage).reshape(-1, 1)
    )

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
def ghi(rd1_file_name, rd2_file_name):
    # TODO: Settle
    BASE_FILE_NAME = "random_seq"
    EXP_CNT = 16
    K_MER_CNT_REP = [EXP_CNT]
    OUTER_DISTANCE = 500
    FRAGMENT_LEN = OUTER_DISTANCE

    # Count k-mers in the random sequence, doubled to represent
    # circular DNA
    start_time = time.time()
    print("Counting k_mers in initial sequence ...", end=" ", flush=True)
    k_mers_seq = count_k_mers_in_seq(initial_seq + initial_seq)
    print("done in {0} s".format(time.time() - start_time), flush=True)

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
    k_mer_cnt_seq = collect_k_mer_cnt(k_mers_seq, seq_cnt=2)
    k_mer_cnt_rd1 = collect_k_mer_cnt(k_mers_rd1, seq_cnt=1)
    k_mer_cnt_rd2 = collect_k_mer_cnt(k_mers_rd2, seq_cnt=1)
    coverage_rd1 = int(np.sum(k_mer_cnt_rd1) / np.sum(k_mer_cnt_seq))
    coverage_rd2 = int(np.sum(k_mer_cnt_rd2) / np.sum(k_mer_cnt_seq))

    # Separate reads into those containing a k-mer with the
    # expected count, and all others
    print("Writing reads based on read count ...", end=" ", flush=True)
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
    trusted_contigs_fnm = os.path.join(
        ssake_out_base_name, ssake_out_base_name + "_scaffolds.fa"
    )
    if os.path.exists(spades_wc_out_dir):
        shutil.rmtree(spades_wc_out_dir)
        command = ru.spades(
            rd1_file_name,
            rd2_file_name,
            spades_wc_out_dir,
            trusted_contigs_fnm=trusted_contigs_fnm,
        )
    print("done in {0} s".format(time.time() - start_time), flush=True)


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
    # Read the multi-FASTA SPAdes output file, and align
    start_time = time.time()
    case = os.path.basename(working_dir)
    print("Aligning {0} assemblies ...".format(case), end=" ", flush=True)
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
