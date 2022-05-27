from pathlib import Path
from time import sleep

import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

import AssemblyUtilities as au
import RunUtilities as ru

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

unq_rps_data = rps_data[
    ~rps_data.duplicated(["Plasmid ID", "Sequence ID"])
].reset_index()
unq_qcs_data = qcs_data.merge(unq_rps_data, on=["Plasmid ID", "Sequence ID"])
unq_qcs_data = unq_qcs_data[
    ~unq_qcs_data.duplicated(["Plasmid ID", "Sequence ID"])
].reset_index()

force = True

for index in range(unq_qcs_data.shape[0]):

    plate = unq_qcs_data.Plate[index]
    well = unq_qcs_data.Well[index]

    working_dir = Path(os.getcwd()) / "assembly" / plate / (well + "/")
    os.makedirs(working_dir, exist_ok=True)

    with au.timing(f"Copying adapters and read files to {working_dir}"):
        adapter_fnm = au.copy_adapters(str(SEQUENCING_BASE_DIR), working_dir)
        # TODO: Add copy only if not exists or force?
        read_fnms = au.copy_reads(str(SEQUENCING_DATA_DIR), plate, well, working_dir)

    for read_fnm in read_fnms:
        kmc_db_fnm = "kmc_" + read_fnm.replace(".fastq.gz", "")
        kmc_dump_fnm = kmc_db_fnm + ".tsv"

        with au.pushd(working_dir):
            if not Path(kmc_dump_fnm).exists() or force:
                with au.timing(f"Counting k_mers in {read_fnm} using kmc"):
                    ru.kmc(
                        [read_fnms[0]],
                        kmc_db_fnm,
                        k_mer_length=25,
                        signature_length=9,
                        count_min=128,  # exclude k-mers with counts less than
                        max_count=1e9,
                        count_max=2048,  # exclude k-mers with counts greater than
                    )
                    ru.kmc_transform(kmc_db_fnm, "dump", kmc_dump_fnm)

            k_mer_cnts = pd.read_csv(kmc_dump_fnm, sep="\t", names=["k_mers", "cnts"])[
                "cnts"
            ].to_numpy()
            rpt_cnt = unq_qcs_data.iloc[index]["Repeat Count"]

            opt_kmeans, opt_s_score, opt_n_clusters = au.find_seq_rcds_for_cnt(
                k_mer_cnts, min_n_clusters=2, max_n_clusters=8
            )

            fig, axs = plt.subplots()
            axs.hist(k_mer_cnts, bins=100)
            axs.set_title(f"k-mer counts distribution\n{well} {plate} {rpt_cnt} {opt_n_clusters}")
            axs.set_xlabel("k-mer counts")
            axs.set_ylabel("frequency")
            plt.show()
            # plt.show(block=False)
            # plt.pause(1)
            # plt.close()
