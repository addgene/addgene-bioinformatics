from pathlib import Path

import os
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

# Consider a sequence
plate = unq_qcs_data.Plate[0]
well = unq_qcs_data.Well[0]

# Create a working directory
working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")
os.makedirs(working_dir, exist_ok=True)

# Copy in adapters and reads
adapter_fnm = au.copy_adapters(str(SEQUENCING_BASE_DIR), working_dir)
read_fnms = au.copy_reads(str(SEQUENCING_DATA_DIR), plate, well, working_dir)

# Process the sequence
with au.pushd(working_dir):

    database_bnms = []
    filtered_fnms = []
    for read_fnm in read_fnms:

        with au.timing(f"Counting k_mers in {read_fnm} using kmc"):
            database_bnm = read_fnm.replace(".fastq.gz", "")
            database_bnms.append(database_bnm)
            ru.kmc(
                [read_fnm],
                database_bnm,
                k_mer_length=25,
                signature_length=9,
                count_min=128,
                max_count=2048,
                count_max=1e9,
            )
            ru.kmc_transform(
                database_bnm, "dump", database_bnm + ".txt", is_sorted=True
            )

        with au.timing(f"Filtering k_mers in {read_fnm} using kmc"):
            filtered_fnm = read_fnm.replace(".fastq.gz", "_filtered.fastq")
            filtered_fnms.append(filtered_fnm)
            ru.kmc_filter(
                database_bnm,
                read_fnm,
                filtered_fnm,
                trim_reads=False,
                hard_mask=False,
                inp_db_count_min=2,
                inp_db_count_max=1e9,
                inp_rd_count_min=64,
                inp_rd_count_max=1e9,
            )

    with au.timing(f"Matching filtered reads"):

        id_sets = []
        for filtered_fnm in filtered_fnms:
            id_set = set()
            with open(filtered_fnm, "r") as f:
                for line in f:
                    if line[0] == "@":
                        id_set.add(line.split(sep=" ")[0])
            id_sets.append(id_set)
        matched_id_set = id_sets[0].intersection(id_sets[1])

        matched_fnms = []
        for filtered_fnm in filtered_fnms:
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
