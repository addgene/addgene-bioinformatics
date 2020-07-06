from argparse import ArgumentParser
from configparser import ConfigParser
import os
import pickle
import time

import test_spades_assembly
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
MAX_ATTEMPTS = 16


if __name__ == "__main__":

    # Add and parse arguments
    parser = ArgumentParser()
    base_dir = os.path.join(os.sep, *os.path.abspath(__file__).split(os.sep)[0:-5])
    parser.add_argument(
        "-s",
        "--sequencing-data-dir",
        default=os.path.join(base_dir, "addgene-sequencing-data", "2018", "FASTQ"),
        help="directory containing sequencing data files",
    )
    parser.add_argument(
        "-a",
        "--assembler-data-dir",
        default=os.path.join(base_dir, "addgene-assembler-data", "results-2020-01-16"),
        help="directory containing assembler data files",
    )
    parser.add_argument(
        "-p",
        "--assembler-data-file",
        default="spades-2020-01-03-174518-pairwise-rl-01.pickle",
        help="file containing assembler data",
    )
    # parser.add_argument(
    #     '-f', '--force-case', type=int,
    #     default=0,
    #     help="force processing of case")
    parser.add_argument(
        "-i",
        "--use-invalid",
        action="store_true",
        help="use invalid cases (no or unaligned circular sequence)",
    )

    OPTIONS = parser.parse_args()

    # Read and parse aligner configuration file
    home_dir = os.path.join(os.sep, *os.path.abspath(__file__).split(os.sep)[0:-4])
    aligner_config_dir = os.path.join(home_dir, "resources")
    aligner_config_file = "pairwise-rl-01.cfg"
    aligner_config = ConfigParser()
    aligner_config.read(os.path.join(aligner_config_dir, aligner_config_file))
    aligner_config["aligner"]["file"] = aligner_config_file

    # Create a pairwise aligner with the specified configuration
    aligner = utilities.create_aligner(aligner_config["aligner"])

    # Load assembly results
    with open(
        os.path.join(OPTIONS.assembler_data_dir, OPTIONS.assembler_data_file), "rb"
    ) as f:
        results = pickle.load(f)

    # Create results directory, if needed
    if not os.path.exists(RESULTS_DIR):
        os.mkdir(RESULTS_DIR)

    # Print header to a file, and stdout
    if OPTIONS.use_invalid:
        output_file = open(__file__.replace(".py", "-invalid.csv"), "w", buffering=1)
    else:
        output_file = open(__file__.replace(".py", "-valid.csv"), "w", buffering=1)
    header = ""
    header += "{:>12s},".format("plate")
    header += "{:>12s},".format("well")
    header += "{:>12s},".format("seq_len")
    header += "{:>12s},".format("exp_spd_scr")
    header += "{:>12s},".format("exp_apc_scr")
    header += "{:>12s},".format("sim_spd_scr")
    header += "{:>12s},".format("sim_apc_scr")
    header += "{:>14s}".format("num_attempts")
    output_file.write(header + "\n")

    # Load alignments, if they exist
    alignments = {}
    pickle_fnm = __file__.replace(".py", ".pickle")
    if os.path.exists(pickle_fnm):
        with open(pickle_fnm, "rb") as pickle_file:
            alignments = pickle.load(pickle_file)

    # Consider each plate and well
    plates = list(results.keys())
    plates.remove("aligner_config")
    wells = list(results[plates[0]].keys())
    for plate in [plates[1]]:
        for well in wells:

            # Skip result if well not found in plate
            if well not in results[plate]:
                continue
            result = results[plate][well]

            # Skip result if no valid QC sequence found
            if "qc" not in result:
                continue
            seq = result["qc"]["sequence"]
            seq_len = len(seq)

            # Skip result if QC sequence has zero length
            if seq_len == 0:
                continue

            # Identify valid cases, those for which an aligned
            # circular sequence exists
            is_valid = False
            exp_spd_scr = 0.0
            exp_apc_scr = 0.0
            if "apc" in result:
                exp_spd_scr = result["spades"]["sequence_score"]
                exp_apc_scr = result["apc"]["sequence_score"]
                if exp_apc_scr == seq_len:
                    is_valid = True

            # Skip cases not selected
            if OPTIONS.use_invalid and is_valid:
                # Using invalid cases, skip valid cases
                continue

            elif not OPTIONS.use_invalid and not is_valid:
                # Using valid cases, skip invalid cases
                continue

            # Skip result if alignment processing complete
            if "results" not in alignments:
                alignments["results"] = []
            if (
                result in alignments["results"]
                and plate in alignments
                and well in alignments[plate]
            ):
                sim_apc_scr = alignments[plate][well]["sim_apc_scr"]
                num_attempts = alignments[plate][well]["num_attempts"]
                if (
                    abs(sim_apc_scr - seq_len) / seq_len <= MAX_DIFFERENCE
                    or num_attempts == MAX_ATTEMPTS
                ):

                    # Format alignments for output
                    seq_len = alignments[plate][well]["seq_len"]
                    exp_spd_scr = alignments[plate][well]["exp_spd_scr"]
                    exp_apc_scr = alignments[plate][well]["exp_apc_scr"]
                    sim_spd_scr = alignments[plate][well]["sim_spd_scr"]
                    sim_apc_scr = alignments[plate][well]["sim_apc_scr"]
                    num_attempts = alignments[plate][well]["num_attempts"]
                    data = ""
                    data += "{:>12s},".format(plate)
                    data += "{:>12s},".format(well)
                    data += "{:>12d},".format(seq_len)
                    data += "{:>12.1f},".format(exp_spd_scr)
                    data += "{:>12.1f},".format(exp_apc_scr)
                    data += "{:>12.1f},".format(sim_spd_scr)
                    data += "{:>12.1f},".format(sim_apc_scr)
                    data += "{:>14d}".format(num_attempts)
                    output_file.write(data + "\n")

                    # If we get here, note processing done
                    print("> Plate {0} well {1} done".format(plate, well))
                    continue

            # If we get here, note processing started
            print("> Plate {0} well {1}:".format(plate, well))

            # Make working directory for well of plate
            working_dir = os.path.join(RESULTS_DIR, plate + "_" + well)
            os.makedirs(working_dir, exist_ok=True)

            # Make directory for assembling simulated reads case
            case_dir = os.path.join(working_dir, "wgsim")
            os.makedirs(case_dir, exist_ok=True)

            # Read data file, if present
            trial_pth = os.path.join(case_dir, "trial.txt")
            if os.path.exists(trial_pth):
                with open(trial_pth, "r") as f:
                    lines = f.readline().replace("\n", "")
                    fields = lines.split(",")
                    sim_spd_scr = float(fields[5])
                    sim_apc_scr = float(fields[6])
                    num_attempts = int(fields[7])
            else:
                sim_spd_scr = 0.0
                sim_apc_scr = 0.0
                num_attempts = 0

            # Attempt to draw simulated reads which assemble correctly
            while (
                abs(sim_apc_scr - seq_len) / seq_len > MAX_DIFFERENCE
                and num_attempts < MAX_ATTEMPTS
            ):
                num_attempts += 1
                print("Trial {0} of {1}".format(num_attempts, MAX_ATTEMPTS))

                # Simulate paired reads of the QC sequence using wgsim,
                # assemble using SPAdes (and apc), and align
                with test_spades_assembly.pushd(case_dir):
                    start_time = time.time()
                    print(
                        "Simulating paired reads from QC sequence using wgsim ...",
                        end=" ",
                        flush=True,
                    )
                    rd1_fnm, rd2_fnm = utilities.simulate_paired_reads_wgsim(
                        seq + seq,
                        WGSIM_BASE_FNM,
                        number_pairs=NUMBER_PAIRS,
                        len_first_read=LEN_FIRST_READ,
                        len_second_read=LEN_SECOND_READ,
                        outer_distance=OUTER_DISTANCE,
                    )
                    print("done in {0} s".format(time.time() - start_time), flush=True)
                try:
                    test_spades_assembly.assemble_using_spades(
                        case_dir, rd1_fnm, rd2_fnm, force=True
                    )
                    (
                        sim_spd_scr,
                        sim_apc_scr,
                    ) = test_spades_assembly.align_assembly_output(
                        aligner, case_dir, seq
                    )
                except Exception as e:
                    print("failed: {0}".format(str(e)), flush=True)

            # Format alignments for output
            data = ""
            data += "{:>12s},".format(plate)
            data += "{:>12s},".format(well)
            data += "{:>12d},".format(seq_len)
            data += "{:>12.1f},".format(exp_spd_scr)
            data += "{:>12.1f},".format(exp_apc_scr)
            data += "{:>12.1f},".format(sim_spd_scr)
            data += "{:>12.1f},".format(sim_apc_scr)
            data += "{:>14d}".format(num_attempts)
            output_file.write(data + "\n")

            # Write data file
            with open(trial_pth, "w") as f:
                f.write(data + "\n")

            # Collect results and alignments, and dump for restarting
            alignments["results"].append(result)
            if plate not in alignments:
                alignments[plate] = {}
            if well not in alignments[plate]:
                alignments[plate][well] = {}
            alignments[plate][well]["seq"] = seq
            alignments[plate][well]["seq_len"] = seq_len
            alignments[plate][well]["exp_spd_scr"] = exp_spd_scr
            alignments[plate][well]["exp_apc_scr"] = exp_apc_scr
            alignments[plate][well]["sim_spd_scr"] = sim_spd_scr
            alignments[plate][well]["sim_apc_scr"] = sim_apc_scr
            alignments[plate][well]["num_attempts"] = num_attempts
            with open(pickle_fnm, "wb") as pickle_file:
                pickle.dump(alignments, pickle_file)
            print("Done")

            # break

        # break
