import itertools
import subprocess
import configparser
from pathlib import Path

OPTIONS = {
    "spades": {"--cov-cutoff": ["auto", "off", "100"], "--careful": [True, False]},
    "unicycler": {
        "--depth-filter": ["0.1", False],
        "‑m": ["normal", "bold", "conservative"],
    },
}
PLATE = "A11967B_sW0154"

for assembler in OPTIONS.keys():
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
                "python",
                "../src/python/toil/PlateAssemblyJob.py",
                f"-o {output_path / PLATE}",
                f"-a {assembler}",
                f"-c {output_path / (assembler + '.ini')}",
                f"-p {PLATE}",
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
        # subprocess.run(command, shell=True)

        # increase the counter
        i += 1
