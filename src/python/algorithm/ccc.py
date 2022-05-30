import logging
import os
from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

import AssemblyUtilities as au
import RunUtilities as ru


SEQUENCING_BASE_DIR = Path("/Volumes/Seagate/addgene-sequencing-data")
SEQUENCING_DATA_DIR = SEQUENCING_BASE_DIR / "2018" / "FASTQ"

QCS_SEQUENCES_FILE = "2018_QC_Sequences.csv"
RPS_SEQUENCES_FILE = (
    "plasmids-sequenced-2018-multiple-features-2022-04-27-ray-will-repeat-features.csv"
)

# Logging configuration
root = logging.getLogger()
if not root.handlers:
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    root.addHandler(ch)

logger = logging.getLogger("ccc")
logger.setLevel(logging.INFO)


# Read QC and repeat sequence data
qcs_data = pd.read_csv(SEQUENCING_BASE_DIR / QCS_SEQUENCES_FILE)
rps_data = pd.read_csv(SEQUENCING_BASE_DIR / RPS_SEQUENCES_FILE)

# Identify the unique intersection of the QC and repeat sequence data
unq_rps_data = rps_data[
    ~rps_data.duplicated(["Plasmid ID", "Sequence ID"])
].reset_index()
unq_qcs_data = qcs_data.merge(unq_rps_data, on=["Plasmid ID", "Sequence ID"])
unq_qcs_data = unq_qcs_data[
    ~unq_qcs_data.duplicated(["Plasmid ID", "Sequence ID"])
].reset_index()

# Consider each sequence
force = False
for i_seq in range(unq_qcs_data.shape[0]):
    plate = unq_qcs_data.iloc[i_seq]["Plate"]
    well = unq_qcs_data.iloc[i_seq]["Well"]

    # Create a working directory
    working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")
    os.makedirs(working_dir, exist_ok=True)
    plot_dir = Path(os.getcwd()) / "plots"
    os.makedirs(plot_dir, exist_ok=True)

    # Copy in adapters and reads
    adapter_fnm = au.copy_adapters(str(SEQUENCING_BASE_DIR), working_dir, force=force)
    raw_read_fnms = au.copy_reads(
        str(SEQUENCING_DATA_DIR), plate, well, working_dir, force=force
    )

    # Preprocess reads
    with au.timing("Preprocessing reads"):
        (
            rd1_trimmed_fnm,
            rd2_trimmed_fnm,
            rd1_normalized_fnm,
            rd2_normalized_fnm,
            rd1_unmerged_fnm,
            rd2_unmerged_fnm,
            rds_merged_fnm,
        ) = au.preprocess_reads_using_bbtools(
            working_dir, raw_read_fnms[0], raw_read_fnms[1], force=force
        )

    # Assemble reads
    with au.timing("Assembling reads"):
        au.assemble_using_spades(
            working_dir,
            raw_read_fnms[0],
            raw_read_fnms[1],
            preprocess=True,
            force=force,
        )

    # Assemble read sets
    count_min = 2
    pp_read_fnms = [
        rd1_unmerged_fnm,
        rd2_unmerged_fnm,
        rds_merged_fnm,
    ]
    with au.pushd(working_dir):
        database_bnms = []
        filtered_fnms = []
        n_clusters = 0
        min_n_clusters = 2
        max_n_clusters = 4
        for read_fnm in pp_read_fnms:

            with au.timing(f"Counting k_mers in {read_fnm} using kmc"):
                database_bnm = read_fnm.replace(".fastq.gz", "")
                database_bnms.append(database_bnm)
                if (
                    not os.path.exists(database_bnm + ".kmc_pre")
                    or not os.path.exists(database_bnm + ".kmc_suf")
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
                if not os.path.exists(dump_fnm) or force:
                    ru.kmc_transform(database_bnm, "dump", dump_fnm, is_sorted=True)

            with au.timing(
                f"Finding optimum number of clusters given counts in {dump_fnm}"
            ):
                if n_clusters != 0:
                    min_n_clusters = n_clusters
                    max_n_clusters = n_clusters
                (
                    k_mer_cnt_min,
                    k_mer_cnt_max,
                ) = au.find_optimum_n_clusters(
                    working_dir, dump_fnm, min_n_clusters=min_n_clusters, max_n_clusters=max_n_clusters, do_plot=False, force=force
                )
                if n_clusters == 0:
                    n_clusters = len(k_mer_cnt_min)

            with au.timing(f"Filtering k_mers in {read_fnm} using kmc"):
                for label in range(n_clusters):
                    filtered_fnm = read_fnm.replace(".fastq.gz", f"_filtered_{label}.fastq")
                    filtered_fnms.append(filtered_fnm)
                    if not os.path.exists(filtered_fnm) or force:
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

        for label in range(n_clusters):
            id_sets = []
            for read_fnm in raw_read_fnms:
                filtered_fnm = read_fnm.replace(".fastq.gz", f"_unmerged_filtered_{label}.fastq")
                id_set = set()
                with open(filtered_fnm, "r") as f:
                    for line in f:
                        if line[0] == "@":
                            id_set.add(line.split(sep=" ")[0])
                id_sets.append(id_set)
            matched_id_set = id_sets[0].intersection(id_sets[1])

            matched_fnms = []
            for read_fnm in raw_read_fnms:
                filtered_fnm = read_fnm.replace(".fastq.gz", f"_unmerged_filtered_{label}.fastq")
                matched_fnm = filtered_fnm.replace("filtered", "matched")
                matched_fnms.append(matched_fnm)
                with open(matched_fnm, "w") as m:
                    with open(filtered_fnm, "r") as f:
                        for line in f:
                            if line[0] == "@":
                                if line.split(sep=" ")[0] in matched_id_set:
                                    do_print = True
                                else:
                                    do_print = False
                            if do_print:
                                m.write(line)

            # Assemble read set
            with au.timing(f"Assembling read set {label} using SPAdes"):
                try:
                    spades_output_dir = au.SPADES_OUTPUT_DIR + f"_{label}"
                    if not os.path.exists(spades_output_dir) or force:
                        ru.spades(
                            rd1_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            rd2_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            spades_output_dir,
                            merged=rds_merged_fnm.replace(".fastq.gz", f"_filtered_{label}.fastq"),
                            cov_cutoff=100,
                        )
                except Exception as ex:
                    logger.error(f"{ex}")

            with au.timing("Assembling paired reads with repeats using SSAKE"):
                try:
                    ssake_output_dir = au.SSAKE_OUTPUT_DIR + f"_{label}"
                    if not os.path.exists(ssake_output_dir) or force:
                        ssake_out_base_name = ssake_output_dir
                        ru.ssake(
                            rd1_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            rd2_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            ssake_out_base_name,
                            au.FRAGMENT_LEN,
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
                except Exception as ex:
                    logger.error(f"{ex}")
