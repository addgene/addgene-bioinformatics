from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class BBDukJob(Job):
    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        config_file_id,
        config_file_name,
        adapters_file_id,
        adapters_file_name,
        parent_rv={},
        *args,
        **kwargs,
    ):

        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short right paired reads
        config_file_id : toil.fileStore.FileID
            id of the file in the file store containing bbduk args
        config_file_name : str
            name of the file in the file store containing bbduk args
        adapters_file_id : toil.fileStore.FileID
            id of the file in the file store containing adapters content
        adapters_file_name : str
            name of the file in the file store containing adapters content
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(BBDukJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name
        self.adapters_file_id = adapters_file_id
        self.adapters_file_name = adapters_file_name
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of bbduk output files
        """
        # Expected output file names
        out1_file_name = "output1.fastq"
        out2_file_name = "output2.fastq"

        try:
            # Read the config file from the file store into the local
            # temporary directory, and parse
            config_file_path = utilities.readGlobalFile(
                fileStore, self.config_file_id, self.config_file_name
            )
            common_config, bbduk_params = utilities.parseConfigFile(
                config_file_path, "bbduk"
            )

            # Read the adapters file from the file store into the local
            # temporary directory
            adapters_file_path = utilities.readGlobalFile(
                fileStore, self.adapters_file_id, self.adapters_file_name
            )

            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, common_config["read_one_file_name"]
            )
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, common_config["read_two_file_name"]
            )

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call bbduk.sh
            # TODO: Specify the container on construction
            image = "ralatsdio/bbtools:v38.90"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))

            # Define BBDuk command
            parameters = [
                "bbduk.sh",
                f"in={read_one_file_path}",
                f"in2={read_two_file_path}",
                f"ref={adapters_file_path}",
                f"out={out1_file_name}",
                f"out2={out2_file_name}",
            ]

            if len(bbduk_params) > 0:
                parameters.extend(bbduk_params)

            logger.info("Using parameters {0}".format(str(parameters)))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
                working_dir=working_dir,
                parameters=parameters,
            )

            # Write the bbduk output files
            # from the local temporary directory into the file store
            out1_file_id = utilities.writeGlobalFile(fileStore, out1_file_name)
            out2_file_id = utilities.writeGlobalFile(fileStore, out2_file_name)

        except Exception as exc:
            # Ensure expected return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            out1_file_id = None
            out2_file_id = None

        # Return file ids and names for export
        bbduk_rv = {
            "bbduk_rv": {
                "out1_file": {
                    "id": out1_file_id,
                    "name": out1_file_name,
                },
                "out2_file": {
                    "id": out2_file_id,
                    "name": out2_file_name,
                },
            }
        }
        bbduk_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(bbduk_rv))
        return bbduk_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification,
    # configuration path and file, adapters path, and output directory, making the
    # output directory if needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
    cmps.extend(["dat", "miscellaneous"])
    parser.add_argument(
        "-d",
        "--data-path",
        default=os.sep + os.path.join(*cmps),
        help="path containing plate and well FASTQ source",
    )
    parser.add_argument(
        "-s", "--source-scheme", default="file", help="scheme used for the source URL"
    )
    parser.add_argument(
        "-l", "--plate-spec", default="A11967A_sW0154", help="the plate specification"
    )
    parser.add_argument(
        "-w", "--well-spec", default="B01", help="the well specification"
    )
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-1]
    parser.add_argument(
        "-c",
        "--config-path",
        default=os.sep + os.path.join(*cmps),
        help="path to a .ini file with args to be passed to the assembler",
    )
    parser.add_argument(
        "-f",
        "--config-file",
        default="Assembler.ini",
        help="path to a .ini file with args to be passed to the assembler",
    )
    parser.add_argument(
        "-ap",
        "--adapters-path",
        default=os.sep + os.path.join(*cmps),
        help="path to a .fa file with adapters to be passed to bbtools",
    )
    parser.add_argument(
        "-a",
        "--adapters-file",
        default="adapters.fa",
        help="path to a .fa file with adapters to be passed to bbtools",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        default=None,
        help="the directory containing all output files",
    )
    options = parser.parse_args()
    if options.output_directory is None:
        options.output_directory = options.plate_spec + "_" + options.well_spec
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)

    # Work within the Toil context manager
    with Toil(options) as toil:
        if not toil.options.restart:

            # Import the local read files into the file store
            read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
                toil,
                options.data_path,
                options.plate_spec,
                [options.well_spec],
                options.source_scheme,
            )

            # Import local config file into the file store
            config_file_id = utilities.importConfigFile(
                toil, os.path.join(options.config_path, options.config_file)
            )

            # Import local adapters file into the file store
            adapters_file_id = utilities.importAdaptersFile(
                toil, os.path.join(options.adapters_path, options.adapters_file)
            )

            # Construct and start the BBDuk job
            bbduk_job = BBDukJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                options.config_file,
                adapters_file_id,
                options.adapters_file,
            )
            bbduk_rv = toil.start(bbduk_job)

        else:

            # Restart the BBDuk job
            bbduk_rv = toil.restart(bbduk_job)

        # Export all BBDuk output files from the file store
        utilities.exportFiles(toil, options.output_directory, bbduk_rv["bbduk_rv"])
