from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class BBMergeJob(Job):
    """
    Accepts paired-end Illumina reads for BBMerge
    (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmerge-guide/)
    """

    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        config_file_id,
        config_file_name,
        chained_job=False,
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
            id of the file in the file store containing bbmerge args
        config_file_name : str
            name of the file in the file store containing bbmerge args
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(BBMergeJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name
        self.chained_job = chained_job
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of bbmerge output files
        """
        # Expected output file names
        out_file_name = "merged.fastq"
        outu1_file_name = "unmerged1.fastq"
        outu2_file_name = "unmerged2.fastq"

        try:
            # Read the config file from the file store into the local
            # temporary directory, and parse
            config_file_path = utilities.readGlobalFile(
                fileStore, self.config_file_id, self.config_file_name
            )
            common_config, bbmerge_params = utilities.parseConfigFile(
                config_file_path, "bbmerge"
            )

            if self.chained_job:
                # Get BBNorm config for input path
                common_config, bbnorm_params = utilities.parseConfigFile(
                    config_file_path, "bbnorm"
                )

                # Read the read files from the file store into the local
                # temporary directory
                read_one_file_path = utilities.readGlobalFile(
                    fileStore,
                    self.read_one_file_id,
                    bbnorm_params["read_one_file_name"],
                )
                read_two_file_path = utilities.readGlobalFile(
                    fileStore,
                    self.read_two_file_id,
                    bbnorm_params["read_two_file_name"],
                )
            else:
                # Read the read files from the file store into the local
                # temporary directory
                read_one_file_path = utilities.readGlobalFile(
                    fileStore,
                    self.read_one_file_id,
                    common_config["read_one_file_name"],
                )
                read_two_file_path = utilities.readGlobalFile(
                    fileStore,
                    self.read_two_file_id,
                    common_config["read_two_file_name"],
                )

            # Read the output filenames from the config file
            outu1_file_name = bbmerge_params["read_one_file_name"]
            outu2_file_name = bbmerge_params["read_two_file_name"]
            out_file_name = bbmerge_params["merged_read_file_name"]

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call bbmerge.sh
            # TODO: Specify the container on construction
            image = "ralatsdio/bbtools:v38.90"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))

            # Define BBMerge command
            parameters = [
                "bbmerge.sh",
                f"in={read_one_file_path}",
                f"in2={read_two_file_path}",
                f"out={out_file_name}",
                f"outu1={outu1_file_name}",
                f"outu2={outu2_file_name}",
            ]

            if len(bbmerge_params) > 0:
                for arg, value in bbmerge_params.items():
                    if (
                        arg != "read_one_file_name"
                        and arg != "read_two_file_name"
                        and arg != "merged_read_file_name"
                    ):
                        parameters.append(value)

            logger.info("Using parameters {0}".format(str(parameters)))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
                working_dir=working_dir,
                parameters=parameters,
            )

            # Write the bbmerge output files
            # from the local temporary directory into the file store
            out_file_id = utilities.writeGlobalFile(fileStore, out_file_name)
            outu1_file_id = utilities.writeGlobalFile(fileStore, outu1_file_name)
            outu2_file_id = utilities.writeGlobalFile(fileStore, outu2_file_name)

        except Exception as exc:
            # Ensure expected return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            out_file_id = None
            outu1_file_id = None
            outu2_file_id = None

        # Return file ids and names for export
        bbmerge_rv = {
            "bbmerge_rv": {
                "merged_file": {
                    "id": out_file_id,
                    "name": out_file_name,
                },
                "outu1_file": {
                    "id": outu1_file_id,
                    "name": outu1_file_name,
                },
                "outu2_file": {
                    "id": outu2_file_id,
                    "name": outu2_file_name,
                },
            }
        }
        bbmerge_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(bbmerge_rv))
        return bbmerge_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification,
    # configuration path and file, and output directory, making the
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

            # Construct and start the BBMerge job
            bbmerge_job = BBMergeJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                options.config_file,
            )
            bbmerge_rv = toil.start(bbmerge_job)

        else:

            # Restart the BBMerge job
            bbmerge_rv = toil.restart(bbmerge_job)

        # Export all BBMerge output files from the file store
        utilities.exportFiles(toil, options.output_directory, bbmerge_rv["bbmerge_rv"])
