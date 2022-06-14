import contextlib
import logging
import os
from pathlib import Path
import shutil
import subprocess
import time

from Bio import Align, SeqIO
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

import RunUtilities as ru


SPADES_OUTPUT_FNM = "contigs.fasta"
SSAKE_FRAGMENT_LEN = 500

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
    plate_dir = list(Path(sequencing_data_dir).glob(plate + "*"))[0].name

    # Follow the read file naming convention
    rd1_fnm = plate_dir.replace("_FASTQ", "") + "_" + well + "_R1_001.fastq.gz"
    rd2_fnm = rd1_fnm.replace("R1", "R2")

    # Copy each file, if needed
    rd_fnms = (rd1_fnm, rd2_fnm)
    for rd_fnm in rd_fnms:
        rd_file_src = Path(sequencing_data_dir) / plate_dir / rd_fnm
        rd_file_dst = Path(working_dir) / rd_fnm
        if not rd_file_dst.exists() or force:
            shutil.copy(rd_file_src, working_dir)

    return rd_fnms


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
    force : boolean
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
    out_fnm = rd1_fnm.replace("_R1", "_Rs").replace(".fastq.gz", "_output.fastq.gz")
    covstats_fnm = rd1_fnm.replace("_R1", "_Rs").replace(".fastq.gz", "_covstats.txt")
    with pushd(working_dir):
        if not Path(out_fnm).exists() or not Path(covstats_fnm).exists() or force:

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


def copy_adapters(sequencing_data_dir, working_dir, force=False):
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
    adapter_file_src = Path(sequencing_data_dir) / adapter_fnm
    adapter_file_dst = Path(working_dir) / adapter_fnm
    if not adapter_file_dst.exists() or force:
        shutil.copy(adapter_file_src, working_dir)

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
    rd1_trimmed_fnm : str
        File name of reads one after BBDuk
    rd2_trimmed_fnm : str
        File name of reads two after BBDuk
    rd1_normalized_fnm : str
        File name of reads one after BBNorm
    rd2_normalized_fnm : str
        File name of reads two after BBNorm
    rd1_unmerged_fnm : str
        File name of reads one after BBMerge
    rd2_unmerged_fnm : str
        File name of reads two after BBMerge
    rds_merged_fnm : str
        File name of merged reads after BBDuk

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
            not Path(rd1_trimmed_fnm).exists()
            or not Path(rd2_trimmed_fnm).exists()
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
            not Path(rd1_normalized_fnm).exists()
            or not Path(rd2_normalized_fnm).exists()
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
            not Path(rds_merged_fnm).exists()
            or not Path(rd1_unmerged_fnm).exists()
            or not Path(rd2_unmerged_fnm).exists()
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


def assemble_using_spades(
    working_dir,
    spades_output_dir,
    apc_base_fnm,
    rd1_fnm,
    rd2_fnm,
    preprocess=True,
    force=False,
):
    """Assemble reads using SPAdes, then circularize using apc.

    Parameters
    ----------
    working_dir : str
        Directory in which the assembly occurs
    spades_output_dir : str
        Directory in which the output is written
    apc_base_fnm : str
        Base file name for apc output
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    preprocess : boolean
        Flag to preprocess reads using BBDuk, BBNorm, and BBMerge
    force : boolean
        Flag to force assembly if output files exist

    Returns
    -------
    spades_seq :  Bio.Seq.Seq
        Sequence object returned by SPAdes
    apc_seq :  Bio.Seq.Seq
        Sequence object returned by apc
    rd1_unmerged_fnm : str
        File name of reads one after optional preprocessing
    rd2_unmerged_fnm : str
        File name of reads two after optional preprocessing
    rds_merged_fnm : str
        File name of Merged reads, after optional preprocessing

    """
    # Optionally preprocess, if needed
    spades_seq = None
    apc_seq = None
    rd1_unmerged_fnm = None
    rd2_unmerged_fnm = None
    rds_merged_fnm = None
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

    # Assemble only if the output directory does not exist
    spades_output_fnm = spades_output_dir + ".fasta"
    spades_output_pth = Path(working_dir) / spades_output_dir / SPADES_OUTPUT_FNM
    apc_output_fnm = apc_base_fnm + ".1.fa"
    with pushd(working_dir):
        if preprocess:

            # Run SPAdes with merged reads, if needed
            if not spades_output_pth.exists() or force:
                with timing("Assembling preprocessed reads using SPAdes"):
                    ru.spades(
                        rd1_unmerged_fnm,
                        rd2_unmerged_fnm,
                        spades_output_dir,
                        merged=rds_merged_fnm,
                        cov_cutoff=100,
                    )

        else:

            # Run SPAdes without merged reads, if needed
            if not spades_output_pth.exists() or force:
                with timing("Assembling original reads using SPAdes"):
                    ru.spades(
                        rd1_fnm,
                        rd2_fnm,
                        spades_output_dir,
                        cov_cutoff=100,
                    )

        # Copy SPAdes output file, and parse
        if spades_output_pth.exists():
            shutil.copy(spades_output_pth, spades_output_fnm)
            spades_seq = [
                seq_rcd for seq_rcd in SeqIO.parse(spades_output_fnm, "fasta")
            ][0].seq

        # Run apc
        if not Path(apc_output_fnm).exists():
            with timing("Circularize using apc"):
                ru.apc(apc_base_fnm, spades_output_fnm)
        logger.info(f"Results of assembling usng SPAdes in directory: {os.getcwd()}")

        # Parse apc output file
        if Path(apc_output_fnm).exists():
            apc_seq = [seq_rcd for seq_rcd in SeqIO.parse(apc_output_fnm, "fasta")][
                0
            ].seq

    return (
        spades_seq,
        apc_seq,
        rd1_unmerged_fnm,
        rd2_unmerged_fnm,
        rds_merged_fnm,
    )


# TODO: Review
def assemble_using_ssake(
        working_dir,
        spades_output_dir,
        apc_base_fnm,
        ssake_scaffolds_fnm,
        rd1_unmerged_fnm,
        rd2_unmerged_fnm,
        rds_merged_fnm,
        force=False
):
    """Assemble reads using SPAdes, then circularize using apc.

    Parameters
    ----------
    working_dir : str
        Directory in which the assembly occurs
    spades_output_dir : str
        Directory in which the output is written
    apc_base_fnm : str
        Base file name for apc output
    ssake_scaffolds_fnm : str
        Name of file containing all scaffolds produced by SSAKE
    rd1_unmerged_fnm : str
        Name of unmerged read one file
    rd2_unmerged_fnm : str
        Name of unmerged read two file
    rds_merged_fnm : str
        Name of merged read one and two file
    force : boolean
        Flag to force assembly if output files exist

    Returns
    -------
    spades_seq :  Bio.Seq.Seq
        Sequence object returned by SPAdes
    apc_seq :  Bio.Seq.Seq
        Sequence object returned by apc
    rd1_matched_m_fnm : str
        File name of reads one matched from minimum k-mer count cluster
    rd2_matched_m_fnm : str
        File name of reads two matched from minimum k-mer count cluster

    """
    spades_seq = None
    apc_seq = None
    rd1_matched_m = None
    rd2_matched_m = None
    with pushd(working_dir):

        # Filter unmerged read files by k-mer count clusters, if
        # needed
        pp_read_fnms = [
            rd1_unmerged_fnm,
            rd2_unmerged_fnm,
        ]
        database_bnms = []
        dump_fnms = []
        count_min = 2
        n_clusters = 0
        min_n_clusters = 2
        max_n_clusters = 4
        for read_fnm in pp_read_fnms:

            # Count k-mers
            with timing(f"Counting k_mers in {read_fnm} using kmc"):
                database_bnm = read_fnm.replace(".fastq.gz", "")
                database_bnms.append(database_bnm)
                if (
                    not Path(database_bnm + ".kmc_pre").exists()
                    or not Path(database_bnm + ".kmc_suf").exists()
                    or force
                ):
                    ru.kmc(
                        [read_fnm],
                        database_bnm,
                        k_mer_length=25,
                        signature_length=9,
                        count_min=count_min,
                        max_count=256,
                        count_max=1e9,
                    )
                dump_fnm = database_bnm + ".txt"
                dump_fnms.append(dump_fnm)
                if not Path(dump_fnm).exists() or force:
                    ru.kmc_transform(database_bnm, "dump", dump_fnm, is_sorted=True)

            # Find k-mer count clusters
            with timing(
                f"Finding optimum number of clusters given counts in {dump_fnm}"
            ):
                if "_R2_" in read_fnm:
                    # Use the optimum number of clusters determined
                    # for the first read file
                    min_n_clusters = n_clusters
                    max_n_clusters = n_clusters
                (k_mer_cnt_min, k_mer_cnt_max,) = find_optimum_n_clusters(
                    dump_fnm,
                    min_n_clusters=min_n_clusters,
                    max_n_clusters=max_n_clusters,
                    do_plot=False,
                    force=force,
                )
                if "_R1_" in read_fnm:
                    # Set the optimum number of clusters, and index of
                    # the cluster with the lowest count determined for
                    # the first read file
                    n_clusters = len(k_mer_cnt_min)
                    i_cluster = np.argmin(np.array(k_mer_cnt_max))

            # Filter reads in k-mer count clusters
            with timing(f"Filtering {read_fnm} using kmc"):
                for label in range(n_clusters):
                    filtered_fnm = read_fnm.replace(
                        ".fastq.gz", f"_filtered_{label}.fastq"
                    )
                    if not Path(filtered_fnm).exists() or force:
                        ru.kmc_filter(
                            database_bnm,
                            read_fnm,
                            filtered_fnm,
                            trim_reads=False,
                            hard_mask=False,
                            inp_db_count_min=k_mer_cnt_min[label],
                            inp_db_count_max=k_mer_cnt_max[label],
                            inp_rd_count_min=2,
                            inp_rd_count_max=1e9,
                        )

        # Match filtered read files
        for label in range(n_clusters):
            filtered_fnms = []
            matched_fnms = []
            for read_fnm in pp_read_fnms:
                filtered_fnm = read_fnm.replace(".fastq.gz", f"_filtered_{label}.fastq")
                filtered_fnms.append(filtered_fnm)
                matched_fnm = read_fnm.replace(".fastq.gz", f"_matched_{label}.fastq")
                matched_fnms.append(matched_fnm)
            with timing(f"Matching paired read files"):
                match_reads(filtered_fnms, matched_fnms)

        # Assemble read sets other than the minumum k-mer count
        # cluster using SSAKE
        with timing(f"Assembling read set {label} using SSAKE"):
            for label in set(range(n_clusters)) - set([i_cluster]):
                try:
                    ssake_output_bnm = f"ssake_{label}"
                    ssake_output_fnm = ssake_output_bnm + "_scaffolds.fa"
                    ssake_output_pth = Path(ssake_output_bnm) / ssake_output_fnm
                    if not ssake_output_pth.exists() or force:
                        if label == 0 and Path(ssake_scaffolds_fnm).exists():
                            os.remove(ssake_scaffolds_fnm)
                        ru.ssake(
                            rd1_unmerged_fnm.replace(
                                ".fastq.gz", f"_matched_{label}.fastq"
                            ),
                            rd2_unmerged_fnm.replace(
                                ".fastq.gz", f"_matched_{label}.fastq"
                            ),
                            ssake_output_bnm,
                            SSAKE_FRAGMENT_LEN,
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

                        # Copy scaffolds to a common file
                        with open(ssake_output_pth, "r") as ifp:
                            with open(ssake_scaffolds_fnm, "a") as ofp:
                                ofp.write(ifp.read())

                except Exception as ex:
                    logger.error(f"{ex}")

        # Assemble unmerged reads from the first k-mer count cluster,
        # all merged reads, and scaffolds using SPAdes and apc
        spades_output_fnm = spades_output_dir + ".fasta"
        spades_output_pth = Path(working_dir) / spades_output_dir / SPADES_OUTPUT_FNM
        apc_output_fnm = apc_base_fnm + ".1.fa"
        with timing(
            f"Assemble unmerged and merged reads, and scaffolds using SPAdes and apc"
        ):
            try:
                if not Path(spades_output_dir).exists() or force:
                    label = i_cluster
                    rd1_matched_m = rd1_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq")
                    rd2_matched_m = rd2_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq")
                    ru.spades(
                        rd1_matched_m,
                        rd2_matched_m,
                        spades_output_dir,
                        merged=rds_merged_fnm,
                        trusted_contigs=ssake_scaffolds_fnm,
                        cov_cutoff=100,
                    )

                # Copy SPAdes output file, and parse
                if spades_output_pth.exists():
                    shutil.copy(spades_output_pth, spades_output_fnm)
                    spades_seq = [
                        seq_rcd for seq_rcd in SeqIO.parse(spades_output_fnm, "fasta")
                    ][0].seq

                # Run apc
                if not Path(apc_output_fnm).exists():
                    with timing("Circularize using apc"):
                        ru.apc(apc_base_fnm, spades_output_fnm)
                logger.info(f"Results of assembling usng SSAKE in directory: {os.getcwd()}")

                # Parse apc output file
                if Path(apc_output_fnm).exists():
                    apc_seq = [seq_rcd for seq_rcd in SeqIO.parse(apc_output_fnm, "fasta")][
                        0
                    ].seq

            except Exception as ex:
                logger.error(f"{ex}")

    return (
        spades_seq,
        apc_seq,
        rd1_matched_m,
        rd2_matched_m,
    )


def find_optimum_n_clusters(
    dump_fnm,
    min_n_clusters=2,
    max_n_clusters=8,
    force=False,
    do_plot=False,
    n_bins=100,
):
    """Use k-means to cluster k-mer counts, selecting the optimum
    number of clusters based on the silhouette score.

    Parameters
    ----------
    dump_fnm : str
        Name of file dumped by kmc
    min_n_clusters : int
        Minimum number of clusters tested
    max_n_clusters : int
        Maximum number of clusters tested
    force : boolean
        Flag to force preprocessing if output files exist
    do_plot : boolean
        Flag to plot histogram with bins color coded by label
    n_bins : int
        Number of bins over which to compute histogram if plotting

    Returns
    -------
    k_mer_cnt_min : int
        Minimum k-mer count for each cluster
    k_mer_cnt_max : int
        Maximum k-mer count for each cluster

    """
    parquet_fnm = dump_fnm.replace(".txt", ".parquet")
    if not Path(parquet_fnm).exists() or force:

        # Read k-mer counts dumped by kmc
        k_mer_cnt_rds = pd.read_csv(dump_fnm, sep="\t", names=["k_mers", "cnts"])

        # Find optimum number of clusters based on the silhouette
        # score
        opt_s_score = -1
        opt_kmeans = None
        cnts = k_mer_cnt_rds["cnts"].to_numpy().reshape(-1, 1)
        with timing("Selecting number of clusters"):
            for n_clusters in range(min_n_clusters, max_n_clusters + 1):
                with timing(f"Fitting and scoring {n_clusters} clusters"):
                    # TODO: Use random state argument
                    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(cnts)
                    s_score = silhouette_score(cnts, kmeans.labels_, random_state=0)
                    logger.info(
                        f"With {n_clusters} clusters found silhouette score {s_score}"
                    )
                    if s_score > opt_s_score:
                        opt_s_score = s_score
                        opt_kmeans = kmeans
        logger.info(
            f"Optimum number of clusters {pd.unique(opt_kmeans.labels_).size} has silhouette score {opt_s_score}"
        )
        k_mer_cnt_rds["labels"] = opt_kmeans.labels_
        k_mer_cnt_rds.to_parquet(parquet_fnm)
    else:
        k_mer_cnt_rds = pd.read_parquet(parquet_fnm)

    # Compute k-mer count minimum and maximum for each label,
    # optionally plotting
    cnts = k_mer_cnt_rds["cnts"]
    labels = k_mer_cnt_rds["labels"]
    if do_plot:
        fig, axs = plt.subplots()
        _n, _bins, patches = axs.hist(cnts, bins=n_bins)
        # TODO: Add more colors
        colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
    n, bins = np.histogram(cnts, bins=n_bins)
    bins = np.array(bins[:-1])
    i_bins = np.array(range(bins.size))
    k_mer_cnt_min = []
    k_mer_cnt_max = []
    for label in range(pd.unique(labels).size):
        k_mer_cnt_min.append(cnts[labels == label].min())
        k_mer_cnt_max.append(cnts[labels == label].max())
        if do_plot:
            for i_bin in i_bins[
                (k_mer_cnt_min[label] <= bins) & (bins <= k_mer_cnt_max[label])
            ]:
                patches[i_bin].set(color=colors[label])
    if do_plot:
        axs.set_title(f"k-mer counts distribution\n{dump_fnm}")
        axs.set_xlabel("k-mer counts")
        axs.set_ylabel("frequency")
        plt.tight_layout()
        plt.savefig(dump_fnm.replace(".txt", ".png"))
        plt.show()

    return k_mer_cnt_min, k_mer_cnt_max


def match_reads(filtered_fnms, matched_fnms):
    """Read filtered one and two read files, find the intersection
    between the identifiers for each read in each file, and write the
    corresponding matched read one and two files.

    Parameters
    ----------
    filtered_fnms : [str]
        Names of filtered read one and two files
    matched_fnms : [str]
        Names of matched read one and two files

    Returns
    -------
    None

    """
    # Check numnber of file names
    if len(filtered_fnms) != len(matched_fnms) and len(filtered_fnms) != 2:
        raise Exception("Number of filtered and matched file names must equal two")

    # Read each filtered read file to create a set of identifiers for
    # each read file, then find the intersection of the sets
    id_sets = []
    for filtered_fnm in filtered_fnms:
        id_set = set()
        with open(filtered_fnm, "r") as f:
            for line in f:
                if line[0] == "@":
                    id_set.add(line.split(sep=" ")[0])
        id_sets.append(id_set)
    matched_id_set = id_sets[0].intersection(id_sets[1])

    # Write the reads in each read file with an identifier in the
    # matched identifier set to a corresponding matched read file
    n_matched = []
    for matched_fnm, filtered_fnm in zip(matched_fnms, filtered_fnms):
        n_matched.append(0)
        with open(matched_fnm, "w") as m:
            with open(filtered_fnm, "r") as f:
                n_group = 3
                do_print = False
                for line in f:
                    if line[0:2] == "@M":
                        if n_group != 3:
                            raise Exception(f"Group prior to {line} did not contain four lines")
                        n_group = 0
                        if line.split(sep=" ")[0] in matched_id_set:
                            do_print = True
                        else:
                            do_print = False
                    else:
                        n_group += 1
                    if do_print:
                        m.write(line)
                        n_matched[-1] += 1
    if n_matched[0] != n_matched[1]:
        raise Exception("Number of matched lines in each read file must be equal")


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
    case = Path(working_dir).name
    with timing(f"Aligning {case} assemblies"):

        # Read the multi-FASTA SPAdes output file, and align
        # TODO: Fix
        output_pth = Path(working_dir) / SPADES_OUTPUT_DIR / SPADES_OUTPUT_FNM
        spd_scr = -1.0
        if Path(output_pth).exists():
            spd_scr = 0.0
            seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
            if len(seq_rcds) > 0:
                spd_scr = aligner.score(seq + seq, seq_rcds[0].seq)

        # Read the FASTA apc output file, and align
        # TODO: Fix
        output_pth = Path(working_dir) / APC_OUTPUT_FNM
        apc_scr = -1.0
        if output_pth.exists():
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
