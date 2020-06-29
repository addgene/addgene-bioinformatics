from argparse import ArgumentParser
from configparser import ConfigParser
import contextlib
import glob
import os
import pickle
import shutil
import time

from Bio import SeqIO

import utilities


RESULTS_DIR = __file__.replace(".py", "")

NUMBER_PAIRS = 25000
LEN_FIRST_READ = 250
LEN_SECOND_READ = 250
OUTER_DISTANCE = 500

SPADES_OUTPUT_DIR = "spades"
SPADES_OUTPUT_FNM = "contigs.fasta"

APC_OUTPUT_DIR = "apc"
APC_BASE_FNM = "apc"
APC_OUTPUT_FNM = APC_BASE_FNM + ".1.fa"

WGSIM_BASE_FNM = "wgsim"
WGSIM_OUTPUT_FNM = WGSIM_BASE_FNM + ".fasta"


@contextlib.contextmanager
def pushd(new_dir):
    """Change to a new working directory, and back.
    """
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    try:
        yield
    finally:
        os.chdir(previous_dir)


def copy_actual_reads(plate, well, working_dir):
    """Copy actual read files for a given plate and well to a working
    directory, if they don't already exist there.
    """
    # Find the plate directory given that we only know the Addgene,
    # but not the seqWell identifier
    plate_dir = os.path.basename(
        glob.glob(
            os.path.join(
                OPTIONS.sequencing_data_dir,
                plate + "*"
            )
        )[0]
    )

    # Follow the read file naming convention
    rd1_fnm = plate_dir.replace(
        "_FASTQ", "") + "_" + well + "_R1_001.fastq.gz"
    rd2_fnm = rd1_fnm.replace("R1", "R2")

    # Copy each file, if needed
    rd_fnms = (rd1_fnm, rd2_fnm)
    for rd_fnm in rd_fnms:
        rd_file_path = os.path.join(
            OPTIONS.sequencing_data_dir,
            plate_dir,
            rd_fnm
        )
        if not os.path.exists(rd_fnm):
            shutil.copy(
                rd_file_path,
                working_dir
            )

    return rd_fnms


def assemble_using_spades(case_dir, rd1_fnm, rd2_fnm, force=False):
    """Assemble reads using SPAdes, then circularize using apc.
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
            command = utilities.spades(
                rd1_fnm, rd2_fnm, SPADES_OUTPUT_DIR, cov_cutoff=100)
            print("done in {0} s".format(time.time() - start_time),
                  flush=True)
            print("Command: {0}".format(command), flush=True)

            # Copy SPAdes output file
            output_pth = os.path.join(
                SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
            if os.path.exists(output_pth):
                shutil.copy(output_pth, ".")

            # Run apc
            start_time = time.time()
            print("Circularize using apc ...", end=" ", flush=True)
            command = utilities.apc(APC_BASE_FNM, SPADES_OUTPUT_FNM)
            print("done in {0} s".format(time.time() - start_time),
                  flush=True)
            print("Command: {0}".format(command), flush=True)


def align_assembly_output(case_dir, seq):

    # Read the multi-FASTA SPAdes output file, and align
    start_time = time.time()
    case = os.path.basename(case_dir)
    print("Aligning {0} assemblies ...".format(case), end=" ", flush=True)
    output_pth = os.path.join(
        case_dir, SPADES_OUTPUT_DIR, SPADES_OUTPUT_FNM)
    spd_scr = 0
    if os.path.exists(output_pth):
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(
            output_pth, "fasta")]
        if len(seq_rcds) > 0:
            spd_scr = aligner.score(seq + seq, seq_rcds[0].seq)

    # Read the FASTA apc output file, and align
    output_pth = os.path.join(case_dir, APC_OUTPUT_FNM)
    apc_scr = 0
    if os.path.exists(output_pth):
        seq_rcds = [seq_rcd for seq_rcd in SeqIO.parse(
            output_pth, "fasta")]
        if len(seq_rcds) > 0:
            apc_scr = aligner.score(seq + seq, seq_rcds[0].seq)
    print("done in {0} s".format(time.time() - start_time), flush=True)
    return spd_scr, apc_scr


if __name__ == "__main__":

    # Add and parse arguments
    parser = ArgumentParser()
    base_dir = os.path.join(
        os.sep, *os.path.abspath(__file__).split(os.sep)[0:-5])
    parser.add_argument(
        '-s', '--sequencing-data-dir',
        default=os.path.join(
            base_dir,
            "addgene-sequencing-data",
            "2018",
            "FASTQ"),
        help="directory containing sequencing data files")
    parser.add_argument(
        '-a', '--assembler-data-dir',
        default=os.path.join(
            base_dir,
            "addgene-assembler-data",
            "results-2020-01-16"),
        help="directory containing assembler data files")
    parser.add_argument(
        '-p', '--assembler-data-file',
        default="spades-2020-01-03-174518-pairwise-rl-01.pickle",
        help="file containing assembler data")
    # parser.add_argument(
    #     '-f', '--force-case', type=int,
    #     default=0,
    #     help="force processing of case")

    OPTIONS = parser.parse_args()

    # Read and parse aligner configuration file
    home_dir = os.path.join(
        os.sep, *os.path.abspath(__file__).split(os.sep)[0:-4])
    aligner_config_dir = os.path.join(home_dir, "resources")
    aligner_config_file = "pairwise-rl-01.cfg"
    aligner_config = ConfigParser()
    aligner_config.read(os.path.join(aligner_config_dir, aligner_config_file))
    aligner_config['aligner']['file'] = aligner_config_file

    # Create a pairwise aligner with the specified configuration
    aligner = utilities.create_aligner(aligner_config['aligner'])

    # Load assembly results
    with open(
        os.path.join(
            OPTIONS.assembler_data_dir, OPTIONS.assembler_data_file
        ), 'rb'
    ) as f:
        results = pickle.load(f)

    # Create results directory, if needed
    if not os.path.exists(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)

    # Print header to a file, and stdout
    output_file = open(__file__.replace(".py", ".csv"), 'w', buffering=1)
    result = ""
    result += "{:>12s},".format("plate")
    result += "{:>12s},".format("well")
    result += "{:>12s},".format("seq_len")
    result += "{:>12s},".format("exp_seq_scr")
    result += "{:>12s},".format("exp_seq_scr")
    result += "{:>12s},".format("act_spd_scr")
    result += "{:>12s},".format("act_apc_scr")
    result += "{:>12s},".format("sim_spd_scr")
    result += "{:>12s}".format("sim_apc_scr")
    output_file.write(result + '\n')

    # Load assemblies, if the exist, and print, or initialize
    assemblies = []
    completed = {}
    pickle_fnm = __file__.replace(".py", ".pickle")
    if os.path.exists(pickle_fnm):
        with open(pickle_fnm, 'rb') as pickle_file:
            assemblies = pickle.load(pickle_file)
            for assembly in assemblies:
                result = ""
                result += "{:>12s},".format(assembly['plate'])
                result += "{:>12s},".format(assembly['well'])
                result += "{:>12d},".format(assembly['seq_len'])
                result += "{:>12.1f},".format(assembly['exp_spd_scr'])
                result += "{:>12.1f},".format(assembly['exp_apc_scr'])
                result += "{:>12.1f},".format(assembly['act_spd_scr'])
                result += "{:>12.1f},".format(assembly['act_apc_scr'])
                result += "{:>12.1f},".format(assembly['sim_spd_scr'])
                result += "{:>12.1f}".format(assembly['sim_apc_scr'])
                output_file.write(result + '\n')

                # Collect completed plates and wells
                plate = assembly['plate']
                well = assembly['well']
                if plate not in completed:
                    completed[plate] = set()
                completed[plate].add(well)

    # Consider each plate and well
    plates = list(results.keys())
    plates.remove('aligner_config')
    wells = list(results[plates[0]].keys())
    for plate in plates:
        for well in wells:

            # Skip result if well found in assemblies
            if plate in completed and well in completed[plate]:
                print("> Plate {0} well {1} done".format(plate, well))
                continue

            # Skip result if well not found in plate
            if well not in results[plate]:
                continue
            assembly = results[plate][well]
            assembly['plate'] = plate
            assembly['well'] = well
            result = ""
            result += "{:>12s},".format(plate)
            result += "{:>12s},".format(well)

            # Skip result if no valid QC assembly found
            if 'qc' not in assembly:
                continue
            seq = assembly['qc']['sequence']
            seq_len = len(seq)
            assembly['seq_len'] = seq_len
            if seq_len == 0:
                continue
            result += "{:>12d},".format(seq_len)

            # Skip result if no circularized assembly found
            if 'apc' not in assembly:
                continue
            exp_spd_scr = assembly['spades']['sequence_score']
            exp_apc_scr = assembly['apc']['sequence_score']
            result += "{:>12.1f},".format(exp_spd_scr)
            result += "{:>12.1f},".format(exp_apc_scr)
            assembly['exp_spd_scr'] = exp_spd_scr
            assembly['exp_apc_scr'] = exp_apc_scr

            # Skip result if sequence score does not equal sequence length
            if (exp_apc_scr != seq_len):
                continue
            print("> Plate {0} well {1}".format(plate, well))

            # Make working directory for well of plate
            working_dir = os.path.join(RESULTS_DIR, plate + "_" + well)
            os.makedirs(working_dir, exist_ok=True)

            # === Actual Reads ===

            # Make directory for assembling actual reads case
            case_dir = os.path.join(working_dir, "actual")
            os.makedirs(case_dir, exist_ok=True)

            # Copy actual reads into place, assemble using SPAdes (and
            # apc), and align
            rd1_fnm, rd2_fnm = copy_actual_reads(plate, well, case_dir)
            assemble_using_spades(case_dir, rd1_fnm, rd2_fnm)
            act_spd_scr, act_apc_scr = align_assembly_output(case_dir, seq)
            result += "{:>12.1f},".format(act_spd_scr)
            result += "{:>12.1f},".format(act_apc_scr)
            assembly['act_spd_scr'] = act_spd_scr
            assembly['act_apc_scr'] = act_apc_scr

            # === Simulated Reads ===

            # Make directory for assembling simulated reads case
            case_dir = os.path.join(working_dir, "wgsim")
            os.makedirs(case_dir, exist_ok=True)

            # Simulate paired reads of the QC sequence using wgsim,
            # assemble using SPAdes (and apc), and align
            if not os.path.exists(
                    os.path.join(case_dir, WGSIM_OUTPUT_FNM)
            ):
                with pushd(case_dir):
                    start_time = time.time()
                    print(
                        "Simulating paired reads from QC sequence using wgsim ...",
                        end=" ", flush=True)
                    rd1_fnm, rd2_fnm = utilities.simulate_paired_reads_wgsim(
                        seq + seq,
                        WGSIM_BASE_FNM,
                        number_pairs=NUMBER_PAIRS,
                        len_first_read=LEN_FIRST_READ,
                        len_second_read=LEN_SECOND_READ,
                        outer_distance=OUTER_DISTANCE)
                    print("done in {0} s".format(time.time() - start_time), flush=True)

            else:
                rd1_fnm = WGSIM_BASE_FNM + "_rd1.fastq"
                rd2_fnm = WGSIM_BASE_FNM + "_rd2.fastq"
            assemble_using_spades(case_dir, rd1_fnm, rd2_fnm)
            sim_spd_scr, sim_apc_scr = align_assembly_output(case_dir, seq)
            result += "{:>12.1f},".format(sim_spd_scr)
            result += "{:>12.1f}".format(sim_apc_scr)
            assembly['sim_spd_scr'] = sim_spd_scr
            assembly['sim_apc_scr'] = sim_apc_scr
            output_file.write(result + '\n')

            # Append current assembly, and dump assemblies for restart
            assemblies.append(assembly)
            with open(pickle_fnm, 'wb') as pickle_file:
                pickle.dump(assemblies, pickle_file)

            # break

        # break
