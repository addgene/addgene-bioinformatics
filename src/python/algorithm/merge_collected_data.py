from pathlib import Path

import pandas as pd


SEQUENCING_BASE_DIR = Path("/Volumes/Seagate/addgene-sequencing-data")
SEQUENCING_DATA_DIR = SEQUENCING_BASE_DIR / "2018" / "FASTQ"

COLLECTED_DATA_HOME = Path("/Users/raymondleclair/Addgene-Spreadsheets")
COLLECTED_DATA_PATH = COLLECTED_DATA_HOME / "2018-Collected-Data.xlsx"


def get_sheet_parquet(sheet_name, usecols=None, header=0):
    sheet_parquet = COLLECTED_DATA_HOME / COLLECTED_DATA_PATH.name.replace(
        ".xlsx", f"-{sheet_name}.parquet"
    )
    if not sheet_parquet.exists():
        sheet = pd.read_excel(
            COLLECTED_DATA_PATH, sheet_name=sheet_name, usecols=usecols, header=header
        )
        sheet.to_parquet(sheet_parquet)
    else:
        sheet = pd.read_parquet(sheet_parquet)
    return sheet


"""
all_seqs_a = pd.read_excel(COLLECTED_DATA_PATH, sheet_name="All-Seqs", header=0)
all_seqs_with_sub_seqs = pd.read_excel(
    COLLECTED_DATA_PATH, sheet_name="All-Seqs-with-Sub-Seqs", header=0
)
all_seqs_without_repeats = pd.read_excel(
    COLLECTED_DATA_PATH, sheet_name="All-Seqs-WITHOUT-Repeats", header=0
)
all_seqs_with_repeats = pd.read_excel(
    COLLECTED_DATA_PATH, sheet_name="All-Seqs-WITH-Repeats", header=0
)
selected_seqs_with_repeats = pd.read_excel(
    COLLECTED_DATA_PATH,
    sheet_name="Selected-Seqs-WITH-Repeats",
    header=0,
    usecols="A:F",
)
all_plas_without_init_seq = pd.read_excel(
    COLLECTED_DATA_PATH, sheet_name="All-Plas-WITHOUT-Init-Seq", header=0
)
all_plas_without_init_with_fini = pd.read_excel(
    COLLECTED_DATA_PATH, sheet_name="All-Plas-WITHOUT-Init-WITH-Fini", header=0
)
all_qc_seqs = pd.read_excel(COLLECTED_DATA_PATH, sheet_name="All-QC-Seqs", header=0)
"""

all_qc_seqs = get_sheet_parquet("All-QC-Seqs")
all_seqs = get_sheet_parquet("All-Seqs")
all_seqs_with_sub_seqs = get_sheet_parquet("All-Seqs-with-Sub-Seqs")
all_seqs_without_repeats = get_sheet_parquet("All-Seqs-WITHOUT-Repeats")
all_seqs_with_repeats = get_sheet_parquet("All-Seqs-WITH-Repeats")
all_plas_without_init_seq = get_sheet_parquet("All-Plas-WITHOUT-Init-Seq")
selected_seqs_with_repeats = get_sheet_parquet(
    "Selected-Seqs-WITH-Repeats", usecols="A:F"
)

# selected_seqs_with_repeat_sub_seqs
columns_a = [
    "Plasmid ID",
    "Sequence ID",
    "Feature Name",
    "Repeat Count",
    "Feature Length",
    # "Notes",
]

# all_seqs_with_sub_seqs
columns_b = [
    "Plasmid ID",
    "Sequence ID",
    "Feature Name",
    "Feature Range Begin",
    "Feature Range End",
    # "Feature Length",
    "Subsequence",
]

# all_qc_seqs
columns_c = [
    "Plasmid ID",
    "Sequence ID",
    "Plate ",
    "Well",
    # "Sequence Length",
    # "Sequence",
]

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

selected_seqs_with_repeat_sub_seqs.to_excel(
    COLLECTED_DATA_HOME / f"2018-Collected-Data-Selected-Seqs-WITH-Repeat-Sub-Seqs.xlsx"
)

selected_seqs_with_repeat_sub_seqs = get_sheet_parquet(
    "Selected-Seqs-WITH-Repeat-Sub-Seqs",
)
