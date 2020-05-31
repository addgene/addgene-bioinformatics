import itertools
import subprocess
import json
from pathlib import Path

OPTIONS = {
    "spades": {"-c": ["auto", "off", "100"], "--careful": [True, False]},
    # "unicycler": {
    #     "--depth-filter": ["0.1", False],
    #     "‑m": ["normal", "bold", "conservative"],
    # },
}

for assembler in OPTIONS.keys():
    i = 0

    # generate every permutation of CLI args
    # modified from https://stackoverflow.com/a/38722093/8068638
    keys, values = zip(*OPTIONS[assembler].items())
    permutations_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]

    # iterate over the various option permutations for the assembler
    for perm in permutations_dicts:

        # hold the args we'll use to call the assembly job
        args = []

        # add the arguments to the list if the key is truthy
        for key, value in perm.items():
            if type(value) == bool and value:
                args.append(key)
            elif type(value) != bool and value:
                args.append((key + "=" + value))

        # create the output path
        output_path = "results/" + assembler + "/" + str(i) + "/"
        if not Path(output_path).exists():
            Path(output_path).mkdir(parents=True)

        command = " ".join(
            [
                "python",
                "../src/python/toil/" + assembler.capitalize() + "Job.py",
                *args,
                "-o=" + output_path,
                "jobstore" + str(i),
            ]
        )

        # for keeping track while running
        print(command)

        # here's where the execution of the command happens
        subprocess.run(command, shell=True)

        # store the args in a parse-friendly way
        with open(output_path + "params.json", "w+") as f:
            json.dump(perm, f)

        # increase the counter
        i += 1
