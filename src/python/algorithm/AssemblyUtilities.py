import contextlib
import glob
import os
import shutil
import time

from Bio import Align
from Bio import SeqIO

import RunUtilities as ru


SPADES_OUTPUT_DIR = "spades"
SPADES_OUTPUT_FNM = "contigs.fasta"

APC_OUTPUT_DIR = "apc"
APC_BASE_FNM = "apc"
APC_OUTPUT_FNM = APC_BASE_FNM + ".1.fa"


@contextlib.contextmanager
def pushd(new_dir):
    """Change to a new working directory, and back.

    Parameters
    ----------
    new_dir : str
        Path of the target directory

    Returns
    -------
    None

    """
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(previous_dir)


def copy_actual_reads(plate, well, working_dir, sequencing_data_dir):
    """Copy actual read files for a given plate and well to a working
    directory, if they don't already exist there.

    Parameters
    ----------
    plate : str
        Plate identifier
    well : str
        Well identifier
    working_dir : str
        Path to working directory
    sequencing_data_dir : str
        Path to directory containing sequencing data

    Returns
    -------
    rd_fnms : tuple(str)
        Tuple containing read file names

    """
    # Find the plate directory given that we only know the Addgene,
    # but not the seqWell identifier
    plate_dir = os.path.basename(
        glob.glob(os.path.join(sequencing_data_dir, plate + "*"))[0]
    )

    # Follow the read file naming convention
    rd1_fnm = plate_dir.replace("_FASTQ", "") + "_" + well + "_R1_001.fastq.gz"
    rd2_fnm = rd1_fnm.replace("R1", "R2")

    # Copy each file, if needed
    rd_fnms = (rd1_fnm, rd2_fnm)
    for rd_fnm in rd_fnms:
        rd_file_path = os.path.join(sequencing_data_dir, plate_dir, rd_fnm)
        if not os.path.exists(rd_fnm):
            shutil.copy(rd_file_path, working_dir)

    return rd_fnms


def assemble_using_spades(case_dir, rd1_fnm, rd2_fnm, force=False):
    """Assemble reads using SPAdes, then circularize using apc.

    Parameters
    ----------
    case_dir : str
        Directory in which the assembly occurs
    rd1_fnm : str
        Name of read one file
    rd2_fnm : str
        Name of read two file
    force=False : boolean
        Flag to force assembly output directory exists

    Returns
    -------
    None

    """
    # Assemble only if the output directory does not exist
    spades_dir = os.path.join(case_dir, SPADES_OUTPUT_DIR)
    if not os.path.exists(spades_dir) or force:
        with pushd(case_dir):

            # if os.path.exists(SPADES_OUTPUT_DIR):
            #     shutil.rmtree(SPADES_OUTPUT_DIR)

            # Run SPAdes
            start_time = time.time()
            print("Assembling using SPAdes ...", end=" ", flush=True)
            command = ru.spades(rd1_fnm, rd2_fnm, SPADES_OUTPUT_DIR, cov_cutoff=100)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)

            # Copy SPAdes output file
            output_pth = os.path.join(SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
            if os.path.exists(output_pth):
                shutil.copy(output_pth, ".")

            # Run apc
            start_time = time.time()
            print("Circularize using apc ...", end=" ", flush=True)
            command = ru.apc(APC_BASE_FNM, SPADES_OUTPUT_FNM)
            print("done in {0} s".format(time.time() - start_time), flush=True)
            print("Command: {0}".format(command), flush=True)
            print("Directory: {0}".format(os.getcwd()))


def align_assembly_output(aligner, case_dir, seq):
    """Align SPAdes and apc assembly output to a sequence.

    Parameters
    ----------
    aligner : Align.PairwiseAligner
        A Biopython pairwise aligner
    case_dir
        Directory containing output of an assembly case
    seq : Bio.Seq
        A Biopython sequence to which to align

    Returns
    -------
    tuple(float)
        Alignment score for the SPAdes and apc output

    """
    # Read the multi-FASTA SPAdes output file, and align
    start_time = time.time()
    case = os.path.basename(case_dir)
    print("Aligning {0} assemblies ...".format(case), end=" ", flush=True)
    output_pth = os.path.join(case_dir, SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
    spd_scr = -1.0
    if os.path.exists(output_pth):
        spd_scr = 0.0
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
        if len(seq_rcds) > 0:
            spd_scr = aligner.score(seq + seq, seq_rcds[0].seq)

    # Read the FASTA apc output file, and align
    output_pth = os.path.join(case_dir, APC_OUTPUT_FNM)
    apc_scr = -1.0
    if os.path.exists(output_pth):
        apc_scr = 0.0
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(output_pth, "fasta")]
        if len(seq_rcds) > 0:
            apc_scr = aligner.score(seq + seq, seq_rcds[0].seq)
    print("done in {0} s".format(time.time() - start_time), flush=True)
    return spd_scr, apc_scr


def create_aligner(config):
    """Creates a pairwise aligner with the specified configuration.

    Parameters
    ----------
    config : dict
        A pairwise aligner configuration

    Returns
    -------
    aligner : Align.PairwiseAligner
        A Biopython pairwise aligner

    """
    # Create a pairwise aligner with the specified configuration
    aligner = Align.PairwiseAligner()
    aligner.match_score = float(config["match_score"])
    aligner.mismatch_score = float(config["mismatch_score"])
    aligner.target_open_gap_score = float(config["target_open_gap_score"])
    aligner.target_extend_gap_score = float(config["target_extend_gap_score"])
    aligner.target_left_open_gap_score = float(config["target_left_open_gap_score"])
    aligner.target_left_extend_gap_score = float(config["target_left_extend_gap_score"])
    aligner.target_right_open_gap_score = float(config["target_right_open_gap_score"])
    aligner.target_right_extend_gap_score = float(
        config["target_right_extend_gap_score"]
    )
    aligner.query_open_gap_score = float(config["query_open_gap_score"])
    aligner.query_extend_gap_score = float(config["query_extend_gap_score"])
    aligner.query_left_open_gap_score = float(config["query_left_open_gap_score"])
    aligner.query_left_extend_gap_score = float(config["query_left_extend_gap_score"])
    aligner.query_right_open_gap_score = float(config["query_right_open_gap_score"])
    aligner.query_right_extend_gap_score = float(config["query_right_extend_gap_score"])
    aligner.mode = config["mode"]

    return aligner
