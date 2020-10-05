import itertools
import subprocess
import configparser
from pathlib import Path
from sys import executable

OPTIONS = {
    "spades": {
        "--cov-cutoff": ["auto", "off", "50", "100", "250"],
        "--careful": [True, False],
        "--plasmid": [True, False],
        "--isolate": [True, False],
    },
    "unicycler": {"--depth_filter": ["0.1", False], "--mode": ["normal", "bold"]},
}
PLATES = [
    "A11935_sW0148",
    "A11938A_sW0144",
    "A11942A_sW0145",
    "A11948_sW0148",
    "A11953A_sW0150",
    "A11956A_sW0152",
    "A11957_sW0148",
    "A11959A_sW0150",
    "A11960A_sW0154",
    "A11966A_sW0152",
]

for assembler in OPTIONS.keys():
    for plate in PLATES:
        i = 0

        # generate every permutation of CLI args
        # modified from https://stackoverflow.com/a/38722093/8068638
        keys, values = zip(*OPTIONS[assembler].items())
        permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

        # iterate over the various option permutations for the assembler
        for perm in permutations_dicts:

            # hold the args we'll use to call the assembly job
            args = {}

            # add the arguments to the list if the key is truthy
            for key, value in perm.items():
                if type(value) == bool and value:
                    args[key] = "true"
                elif type(value) != bool and value:
                    args[key] = value

            # create the output path
            output_path = Path("results/" + assembler + "/" + str(i) + "/")
            if not output_path.exists():
                output_path.mkdir(parents=True)

            command = " ".join(
                [
                    executable,
                    "../jobs/PlateAssemblyJob.py",
                    f"-o {output_path / plate}",
                    f"-a {assembler}",
                    f"-c {output_path / (assembler + '.ini')}",
                    f"-p {plate}",
                    "-s",
                    "s3",
                    "-d",
                    "addgene-sequencing-data/2018/FASTQ",
                    f"jobstore{str(i)}",
                ]
            )

            # for keeping track while running
            print(command)

            # store the args in a parse-friendly way
            config = configparser.ConfigParser()
            config[assembler] = args
            with open(output_path / f"{assembler}.ini", "w+") as f:
                config.write(f)

            # here's where the execution of the command happens
            subprocess.run(command, shell=True)

            # increase the counter
            i += 1
