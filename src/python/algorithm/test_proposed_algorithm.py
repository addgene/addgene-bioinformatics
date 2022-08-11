from configparser import ConfigParser
import logging
import os
from pathlib import Path

import pandas as pd

import AssemblyUtilities as au


SEQUENCING_BASE_DIR = Path("/Volumes/Seagate/addgene-sequencing-data")
SEQUENCING_DATA_DIR = SEQUENCING_BASE_DIR / "2018" / "FASTQ"

COLLECTED_DATA_HOME = Path("/Users/raymondleclair/Addgene-Spreadsheets")
COLLECTED_DATA_PATH = COLLECTED_DATA_HOME / "2018-Collected-Data.xlsx"


# Logging configuration
root = logging.getLogger()
if not root.handlers:
    ch = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    root.addHandler(ch)

logger = logging.getLogger("test_proposed_algorithm")
logger.setLevel(logging.INFO)


def merge_collected_data():
    """Merge collected data provided by Addgene to produce a dataframe
    with plasmid and sequence identifier, plate and well identifier,
    and repeated subsequences.

    Parameters
    ----------
    None

    Returns
    -------
    selected_seqs_with_repeat_sub_seqs : pd.DataFrame
        Pandas DataFrame containing merged collected data
    """
    # Load all sheets provided by Addgene from the collected data
    # spreadsheet
    logger.info("Loading All-QC-Seqs")
    all_qc_seqs = au.get_sheet_parquet(
        COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-QC-Seqs"
    )
    # all_seqs = au.get_sheet_parquet(
    #     COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Seqs"
    # )
    logger.info("Loading All-Seqs-with-Sub-Seqs")
    all_seqs_with_sub_seqs = au.get_sheet_parquet(
        COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Seqs-with-Sub-Seqs"
    )
    # all_seqs_without_repeats = au.get_sheet_parquet(
    #     COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Seqs-WITHOUT-Repeats"
    # )
    # all_seqs_with_repeats = au.get_sheet_parquet(
    #     COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Seqs-WITH-Repeats"
    # )
    # all_plas_without_init_seq = au.get_sheet_parquet(
    #     COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Plas-WITHOUT-Init-Seq"
    # )
    logger.info("Loading Selected-Seqs-WITH-Repeats")
    selected_seqs_with_repeats = au.get_sheet_parquet(
        COLLECTED_DATA_HOME,
        COLLECTED_DATA_PATH,
        "Selected-Seqs-WITH-Repeats",
        usecols="A:F",
    )

    # Define columns needed for merge in the tab ...
    # .... selected_seqs_with_repeat_sub_seqs
    columns_a = [
        "Plasmid ID",
        "Sequence ID",
        "Feature Name",
        "Repeat Count",
        "Feature Length",
        # "Notes",
    ]
    # ... all_seqs_with_sub_seqs
    columns_b = [
        "Plasmid ID",
        "Sequence ID",
        "Feature Name",
        "Feature Range Begin",
        "Feature Range End",
        # "Feature Length",
        "Subsequence",
    ]
    # ... all_qc_seqs
    columns_c = [
        "Plasmid ID",
        "Sequence ID",
        "Plate",
        "Well",
        # "Sequence Length",
        # "Sequence",
    ]

    # Perform merge
    logger.info("Merging collected data")
    selected_seqs_with_repeat_sub_seqs = (
        selected_seqs_with_repeats[columns_a]
        .merge(
            all_seqs_with_sub_seqs[columns_b],
            how="left",
            on=["Plasmid ID", "Sequence ID", "Feature Name"],
        )
        .merge(
            all_qc_seqs[columns_c],
            how="left",
            on=["Plasmid ID", "Sequence ID"],
        )
    )

    # Save merge result to separate Excel file
    logger.info("Saving merged data")
    selected_seqs_with_repeat_sub_seqs.to_excel(
        COLLECTED_DATA_HOME
        / "2018-Collected-Data-Selected-Seqs-WITH-Repeat-Sub-Seqs.xlsx"
    )

    return all_qc_seqs, selected_seqs_with_repeat_sub_seqs


def main():

    # Read and parse configuration file
    home_dir = Path(__file__).parents[3]
    config_dir = home_dir / "resources"
    config_file = config_dir / "pairwise-rl-01.cfg"
    config = ConfigParser()
    config.read(config_file)
    config["aligner"]["file"] = str(config_file)

    # Read QC and selected sequences with repeat subsequences
    qc_seqs, sel_seqs = merge_collected_data()

    # Create a pairwise aligner with the specified configuration
    aligner = au.create_aligner(dict(config["aligner"]))

    # Consider each unique sequence
    force = False
    unq_seqs = sel_seqs[
        ~sel_seqs.duplicated(["Plasmid ID", "Sequence ID"])
    ].reset_index()
    for i_seq in range(unq_seqs.shape[0]):
        sequence_id = unq_seqs.iloc[i_seq]["Sequence ID"]
        plate = unq_seqs.iloc[i_seq]["Plate"]
        well = unq_seqs.iloc[i_seq]["Well"]

        # Get QC and subsequences
        qc_seq = qc_seqs.loc[qc_seqs["Sequence ID"] == sequence_id]["Sequence"].iloc[0]
        sub_seqs = [
            sub_seq
            for sub_seq in sel_seqs.loc[sel_seqs["Sequence ID"] == sequence_id][
                "Subsequence"
            ]
        ]

        # Create a directory for plot files
        plot_dir = Path(os.getcwd()) / "plots"
        os.makedirs(plot_dir, exist_ok=True)

        # Create a working directory
        working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")
        os.makedirs(working_dir, exist_ok=True)

        # Copy in reads and adapters
        with au.timing(
            f"Copying in reads for plate {plate} and well {well} and adapters"
        ):
            raw_read_fnms = au.copy_reads(
                str(SEQUENCING_DATA_DIR), plate, well, working_dir, force=force
            )
            adapter_fnm = au.copy_adapters(
                SEQUENCING_BASE_DIR, working_dir, force=force
            )

        # Assemble reads using SPAdes
        with au.timing("Assembling reads using SPAdes"):
            (
                spades_seq_current,
                apc_seq_current,
                rd1_unmerged_fnm,
                rd2_unmerged_fnm,
                rds_merged_fnm,
            ) = au.assemble_using_spades(
                working_dir,
                "spades_current",
                "apc_current",
                raw_read_fnms[0],
                raw_read_fnms[1],
                preprocess=True,
                force=force,
            )

        # Align result to QC sequence
        print("Start here")
        spades_scr_current = 0.0
        if spades_seq_current is not None:
            spades_scr_current = aligner.score(qc_seq + qc_seq, spades_seq_current)
        apc_scr_current = 0.0
        if apc_seq_current is not None:
            apc_scr_current = aligner.score(qc_seq + qc_seq, apc_seq_current)

        # Assemble reads using SSAKE
        with au.timing("Assembling reads using SSAKE"):
            (
                spades_seq_proposed,
                apc_seq_proposed,
                _rd1_matched_m,
                _rd2_matched_m,
            ) = au.assemble_using_ssake(
                working_dir,
                "spades_proposed",
                "apc_proposed",
                "ssake_scaffolds.fasta",
                rd1_unmerged_fnm,
                rd2_unmerged_fnm,
                rds_merged_fnm,
                force=force,
            )

        # Align result to QC sequence
        spades_scr_proposed = 0.0
        if spades_seq_proposed is not None:
            spades_scr_proposed = aligner.score(qc_seq + qc_seq, spades_seq_proposed)
        apc_scr_proposed = 0.0
        if apc_seq_proposed is not None:
            apc_scr_proposed = aligner.score(qc_seq + qc_seq, apc_seq_proposed)

        # TODO: Count repeats in SSAKE scaffolds


if __name__ == "__main__":
    main()
