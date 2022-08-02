
import logging
import os
from pathlib import Path

import pandas as pd

import AssemblyUtilities as au


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

    # Assemble reads using SPAdes
    with au.timing("Assembling reads using SPAdes"):
        (
            spades_seq,
            apc_seq,
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

    # TODO: Align result to QC sequence

    # Assemble reads using SSAKE
    with au.timing("Assembling reads using SSAKE"):
        au.assemble_using_ssake(
            working_dir,
            "spades_proposed",
            "apc_proposed",
            "ssake_scaffolds.fasta",
            rd1_unmerged_fnm,
            rd2_unmerged_fnm,
            rds_merged_fnm,
            force=force,
        )

    # TODO: Align result to QC sequence
