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

    # Copy in reads and adapters
    raw_read_fnms = au.copy_reads(
        str(SEQUENCING_DATA_DIR), plate, well, working_dir, force=force
    )
    adapter_fnm = au.copy_adapters(SEQUENCING_BASE_DIR, working_dir, force=force)

    # Preprocess reads
    # with au.timing("Preprocessing reads"):
    #     (
    #         rd1_trimmed_fnm,
    #         rd2_trimmed_fnm,
    #         rd1_normalized_fnm,
    #         rd2_normalized_fnm,
    #         rd1_unmerged_fnm,
    #         rd2_unmerged_fnm,
    #         rds_merged_fnm,
    #     ) = au.preprocess_reads_using_bbtools(
    #         working_dir, raw_read_fnms[0], raw_read_fnms[1], force=force
    #     )

    # Assemble reads with preprocessing
    with au.timing("Assembling reads"):
        (
            spades_seq,
            apc_seq,
            rd1_unmerged_fnm,
            rd2_unmerged_fnm,
            rds_merged_fnm,
        ) = au.assemble_using_spades(
            working_dir,
            "spades_0",
            "apc_0",
            raw_read_fnms[0],
            raw_read_fnms[1],
            preprocess=True,
            force=force,
        )

    # TODO: Align result to QC sequence

    # TODO: Move all of the following to AssemblyUtilities

    with au.pushd(working_dir):

        # Filter unmerged read files by k-mer count clusters
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

            with au.timing(f"Counting k_mers in {read_fnm} using kmc"):
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

            with au.timing(
                f"Finding optimum number of clusters given counts in {dump_fnm}"
            ):
                if "_R2_" in read_fnm:
                    # Use the optimum number of clusters determined for the first read file
                    min_n_clusters = n_clusters
                    max_n_clusters = n_clusters
                (
                    k_mer_cnt_min,
                    k_mer_cnt_max,
                ) = au.find_optimum_n_clusters(
                    dump_fnm, min_n_clusters=min_n_clusters, max_n_clusters=max_n_clusters, do_plot=False, force=force
                )
                if "_R1_" in read_fnm:
                    # Set the optimum number of clusters, and index of the cluster with the lowest count determined for the first read file
                    n_clusters = len(k_mer_cnt_min)
                    i_cluster = np.argmin(np.array(k_mer_cnt_max))

            with au.timing(f"Filtering {read_fnm} using kmc"):
                for label in range(n_clusters):
                    filtered_fnm = read_fnm.replace(".fastq.gz", f"_filtered_{label}.fastq")
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
            with au.timing(f"Matching paired read files"):
                au.match_reads(filtered_fnms, matched_fnms)

        # Assemble read sets using SSAKE
        with au.timing(f"Assembling read set {label} using SSAKE"):
            ssake_scaffolds_fnm = "ssake_scaffolds.fa"
            for label in set(range(n_clusters)) - set([i_cluster]):
                try:
                    ssake_output_bnm = f"ssake_{label}"
                    ssake_output_fnm = ssake_output_bnm + "_scaffolds.fa"
                    ssake_output_pth = Path(ssake_output_bnm) / ssake_output_fnm
                    if not ssake_output_pth.exists() or force:
                        if label == 0 and Path(ssake_scaffolds_fnm).exists():
                            os.remove(ssake_scaffolds_fnm)
                        ru.ssake(
                            rd1_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            rd2_unmerged_fnm.replace(".fastq.gz", f"_matched_{label}.fastq"),
                            ssake_output_bnm,
                            au.SSAKE_FRAGMENT_LEN,
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

        # Assemble unmerged reads from the first k-mer count cluster, all merged reads, and scaffolds using SPAdes and apc
        with au.timing(f"Assemble unmerged and merged reads, and scaffolds using SPAdes and apc"):
            try:
                spades_output_dir = "spades_s"
                if not Path(spades_output_dir).exists() or force:
                    ru.spades(
                        rd1_unmerged_fnm,
                        rd2_unmerged_fnm,
                        spades_output_dir,
                        merged=rds_merged_fnm,
                        trusted_contigs=ssake_scaffolds_fnm,
                        cov_cutoff=100,
                    )
            except Exception as ex:
                logger.error(f"{ex}")

            # TODO: Align result to QC sequence
