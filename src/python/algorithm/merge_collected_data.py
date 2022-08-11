from pathlib import Path

import AssemblyUtilities as au


SEQUENCING_BASE_DIR = Path("/Volumes/Seagate/addgene-sequencing-data")
SEQUENCING_DATA_DIR = SEQUENCING_BASE_DIR / "2018" / "FASTQ"

COLLECTED_DATA_HOME = Path("/Users/raymondleclair/Addgene-Spreadsheets")
COLLECTED_DATA_PATH = COLLECTED_DATA_HOME / "2018-Collected-Data.xlsx"


def main():
    # Load all sheets provided by Addgene from the collected data
    # spreadsheet
    all_qc_seqs = au.get_sheet_parquet(
        COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-QC-Seqs"
    )
    # all_seqs = au.get_sheet_parquet(
    #     COLLECTED_DATA_HOME, COLLECTED_DATA_PATH, "All-Seqs"
    # )
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
        "Plate ",
        "Well",
        # "Sequence Length",
        # "Sequence",
    ]

    # Perform merge
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
    selected_seqs_with_repeat_sub_seqs.to_excel(
        COLLECTED_DATA_HOME
        / "2018-Collected-Data-Selected-Seqs-WITH-Repeat-Sub-Seqs.xlsx"
    )

    # Get merge result from collected data Excel file
    # Note: The following requires a manual copy
    try:
        _selected_seqs_with_repeat_sub_seqs = au.get_sheet_parquet(
            COLLECTED_DATA_HOME,
            COLLECTED_DATA_PATH,
            "Sel-Seqs-WITH-Repeat-Sub-Seqs",
        )
    except Exception as ex:
        print(
            f"Could not get sheet parquet for 'Selected-Seqs-WITH-Repeat-Sub-Seqs': {ex}"
        )

    return selected_seqs_with_repeat_sub_seqs


if __name__ == "__main__":
    sequences = main()
