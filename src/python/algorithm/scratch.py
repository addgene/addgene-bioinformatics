from pathlib import Path

import os
import pandas as pd

import AssemblyUtilities as au


SEQUENCING_BASE_DIR = Path("/Volumes/Seagate/addgene-sequencing-data")
SEQUENCING_DATA_DIR = SEQUENCING_BASE_DIR / "2018" / "FASTQ"

QCS_SEQUENCES_FILE = "2018_QC_Sequences.csv"
RPS_SEQUENCES_FILE = (
    "plasmids-sequenced-2018-multiple-features-2022-04-27-ray-will-all-sequences.csv"
)
RPS_SEQUENCES_FILE = (
    "plasmids-sequenced-2018-multiple-features-2022-04-27-ray-will-mutiple-features.csv"
)
RPS_SEQUENCES_FILE = (
    "plasmids-sequenced-2018-multiple-features-2022-04-27-ray-will-repeat-features.csv"
)

qcs_data = pd.read_csv(SEQUENCING_BASE_DIR / QCS_SEQUENCES_FILE)
rps_data = pd.read_csv(SEQUENCING_BASE_DIR / RPS_SEQUENCES_FILE)

unq_rps_data = rps_data[~rps_data.duplicated(["Plasmid ID", "Sequence ID"])].reset_index()
unq_qcs_data = qcs_data.merge(unq_rps_data, on=["Plasmid ID", "Sequence ID"])
unq_qcs_data = unq_qcs_data[~unq_qcs_data.duplicated(["Plasmid ID", "Sequence ID"])].reset_index()

plate = unq_qcs_data.Plate[0]
well = unq_qcs_data.Well[0]

plate = "A11967A"
well = "B01"

working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")

os.makedirs(working_dir, exist_ok=True)

adapter_fnm = au.copy_adapters(str(SEQUENCING_DATA_DIR), working_dir) 
read_fnms = au.copy_reads(str(SEQUENCING_DATA_DIR), plate, well, working_dir)

au.assemble_using_spades(working_dir, read_fnms[0], read_fnms[1], preprocess=True, force=True)
