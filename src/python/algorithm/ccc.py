import os
from pathlib import Path

from matplotlib import pyplot as plt
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

# Consider each sequence
for index in range(unq_qcs_data.shape[0]):
    plate = unq_qcs_data.Plate[index]
    well = unq_qcs_data.Well[index]

    # Create a working directory
    working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")
    os.makedirs(working_dir, exist_ok=True)
    plot_dir = Path(os.getcwd()) / "plots"
    os.makedirs(plot_dir, exist_ok=True)

    # Copy in adapters and reads
    adapter_fnm = au.copy_adapters(str(SEQUENCING_BASE_DIR), working_dir)
    raw_read_fnms = au.copy_reads(str(SEQUENCING_DATA_DIR), plate, well, working_dir)

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
            working_dir, raw_read_fnms[0], raw_read_fnms[1], force=True
        )
    pp_read_fnms = [
        rd1_unmerged_fnm,
        rd2_unmerged_fnm,
        rds_merged_fnm,
    ]

    count_min = 2

    # Process the sequence
    with au.pushd(working_dir):

        database_bnms = []
        filtered_fnms = []
        for read_fnm in pp_read_fnms:

            with au.timing(f"Counting k_mers in {read_fnm} using kmc"):
                database_bnm = read_fnm.replace(".fastq.gz", "")
                database_bnms.append(database_bnm)
                dump_fnm = database_bnm + ".txt"
                ru.kmc(
                    [read_fnm],
                    database_bnm,
                    k_mer_length=25,
                    signature_length=9,
                    count_min=count_min,
                    max_count=256,
                    count_max=1e9,
                )
                ru.kmc_transform(
                    database_bnm, "dump", dump_fnm, is_sorted=True
                )
                k_mer_cnts = pd.read_csv(dump_fnm, sep="\t", names=["k_mers", "cnts"])[
                    "cnts"
                    ].to_numpy()
                rpt_cnt = unq_qcs_data.iloc[index]["Repeat Count"]
                fig, axs = plt.subplots()
                axs.hist(k_mer_cnts, bins=100)
                axs.set_title(f"k-mer counts distribution\n{read_fnm} {rpt_cnt}")
                axs.set_xlabel("k-mer counts")
                axs.set_ylabel("frequency")
                plt.tight_layout()
                plt.savefig(plot_dir / dump_fnm.replace(".txt", ".png"))
                plt.show(block=False)
                plt.pause(1)
                plt.close()

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
                    inp_rd_count_min=count_min,
                    inp_rd_count_max=1e9,
                )
