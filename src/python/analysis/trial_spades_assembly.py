from argparse import ArgumentParser
from configparser import ConfigParser
import os
import pickle
import time

from test_spades_assembly import align_assembly_output
from test_spades_assembly import assemble_using_spades
from test_spades_assembly import pushd
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

MAX_DIFFERENCE = 0.10
MAX_ATTEMPTS = 8


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
    result += "{:>12s},".format("sim_spd_scr")
    result += "{:>12s},".format("sim_apc_scr")
    result += "{:>14s}".format("num_attempts")
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
                result += "{:>12.1f},".format(assembly['sim_spd_scr'])
                result += "{:>12.1f},".format(assembly['sim_apc_scr'])
                result += "{:>14d}".format(assembly['num_attempts'])
                output_file.write(result + '\n')

                # Collect completed plates and wells
                plate = assembly['plate']
                if plate not in completed:
                    completed[plate] = {}
                well = assembly['well']
                seq_len = assembly['seq_len']
                sim_apc_scr = assembly['sim_apc_scr']
                num_attempts = assembly['num_attempts']
                if (abs(sim_apc_scr - seq_len) / seq_len <= MAX_DIFFERENCE
                    or num_attempts == MAX_ATTEMPTS):
                    completed[plate][well] = True
                else:
                    completed[plate][well] = False

    # Consider each plate and well
    plates = list(results.keys())
    plates.remove('aligner_config')
    wells = list(results[plates[0]].keys())
    for plate in [plates[1]]:
        for well in wells:

            # Skip result if well found in assemblies
            if (plate in completed
                and well in completed[plate]
                and completed[plate][well]):
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
            assembly['exp_spd_scr'] = exp_spd_scr
            assembly['exp_apc_scr'] = exp_apc_scr
            result += "{:>12.1f},".format(exp_spd_scr)
            result += "{:>12.1f},".format(exp_apc_scr)

            # Skip result if sequence score does not equal sequence length
            if (exp_apc_scr != seq_len):
                continue
            print("> Plate {0} well {1}:".format(plate, well))

            # Make working directory for well of plate
            working_dir = os.path.join(RESULTS_DIR, plate + "_" + well)
            os.makedirs(working_dir, exist_ok=True)

            # Make directory for assembling simulated reads case
            case_dir = os.path.join(working_dir, "wgsim")
            os.makedirs(case_dir, exist_ok=True)

            # Attempt to draw simulated reads which assemble correctly
            sim_spd_scr, sim_apc_scr = 0, 0
            num_attempts = 0
            while (abs(sim_apc_scr - seq_len) / seq_len > MAX_DIFFERENCE
                   and num_attempts < MAX_ATTEMPTS):
                num_attempts += 1
                print("Trial {0} of {1}".format(num_attempts, MAX_ATTEMPTS))
                
                # Simulate paired reads of the QC sequence using wgsim,
                # assemble using SPAdes (and apc), and align
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
                    print("done in {0} s".format(time.time() - start_time),
                          flush=True)
                sim_spd_scr, sim_apc_scr = 0, 0
                try:
                    assemble_using_spades(case_dir, rd1_fnm, rd2_fnm, force=True)
                    sim_spd_scr, sim_apc_scr = align_assembly_output(aligner, case_dir, seq)
                except Exception as e:
                    print("failed: {0}".format(str(e)), flush=True)

            # Collect results for restarting, and write to output file
            assembly['sim_spd_scr'] = sim_spd_scr
            assembly['sim_apc_scr'] = sim_apc_scr
            assembly['num_attempts'] = num_attempts
            result += "{:>12.1f},".format(sim_spd_scr)
            result += "{:>12.1f},".format(sim_apc_scr)
            result += "{:>14d}".format(num_attempts)
            output_file.write(result + '\n')

            # Append current assembly, and dump assemblies for restart
            assemblies.append(assembly)
            with open(pickle_fnm, 'wb') as pickle_file:
                pickle.dump(assemblies, pickle_file)

            # break

        # break
