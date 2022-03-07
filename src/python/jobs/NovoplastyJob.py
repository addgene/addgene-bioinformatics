from argparse import ArgumentParser
import logging
import os
import gzip

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class NovoplastyJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using NOVOPlasty.
    """

    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        config_file_id,
        config_file_name,
        chained_job=False,
        seed_file_id=False,
        parent_rv={},
        *args,
        **kwargs
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
            id of the file in the file store containing assembler args
        config_file_name : str
            name of the file in the file store containing assembler args
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(NovoplastyJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name
        self.chained_job = chained_job
        self.parent_rv = parent_rv
        self.seed_file_id = seed_file_id

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of log and corrections files, and
            contigs FASTA file
        """
        # Expected output file names
        project_name = "Toil"
        log_file_name = "log_{0}.txt".format(project_name)
        contigs_file_name = "Circularized_assembly_1_{0}.fasta".format(project_name)

        try:
            # Read the config files from the file store into the local
            # temporary directory, and parse
            config_file_path = utilities.readGlobalFile(
                fileStore, self.config_file_id, self.config_file_name
            )
            common_config, assembler_params = utilities.parseConfigFile(
                config_file_path, "novoplasty"
            )
            common_config, bbnorm_params = utilities.parseConfigFile(
                config_file_path, "bbnorm"
            )

            # Read the read files from the file store into the local
            # temporary directory
            if self.chained_job:
                read_one_file_name = bbnorm_params["read_one_file_name"]
                read_two_file_name = bbnorm_params["read_two_file_name"]
            else:
                read_one_file_name = common_config["read_one_file_name"]
                read_two_file_name = common_config["read_two_file_name"]
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, read_one_file_name
            )
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, read_two_file_name
            )

            # Select the first read sequence as the seed, and write it
            # into the local temporary directory
            working_dir = fileStore.localTempDir

            # Check if seed file exists; I.E., it's been imported
            if self.seed_file_id:
                # Read existing seed file
                seed_file_path = utilities.readGlobalFile(
                    fileStore, self.seed_file_id, os.path.basename(self.seed_file_id)
                )
            else:
                seed_file_path = os.path.join(
                    working_dir, assembler_params["seed_file_name"]
                )

                # Create seed file from first read
                with open(seed_file_path, "w+") as f:
                    with open(read_one_file_path, "rt") as g:
                        do_write = False
                        for line in g:
                            if line[0] == "@":
                                do_write = True
                            if line[0] == "+":
                                break
                            if do_write:
                                f.write(line)

            # Write the NOVOPlasty config file into the local temporary
            # directory
            logger.info(
                "Handling configuration file {0}".format(
                    assembler_params["config_file_name"]
                )
            )
            with open(
                os.path.join(working_dir, assembler_params["config_file_name"]), "w+"
            ) as f:
                config = """Project:
-----------------------
Project name          = {project_name}
Type                  = mito
Genome Range          = 1500-60000
K-mer                 = 121
Max memory            = 3
Extended log          = 0
Save assembled reads  = no
Seed Input            = {seed_file_path}
Extend seed directly  = no
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 251
Insert size           = 500
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = {read_one_file_path}
Reverse reads         = {read_two_file_path}

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
                 """.format(
                    project_name=project_name,
                    seed_file_path=seed_file_path,
                    read_one_file_path=read_one_file_path,
                    read_two_file_path=read_two_file_path,
                )
                f.write(config)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call NOVOPlasty
            # TODO: Specify the container on construction
            image = "ralatsdio/novoplasty:v4.2"
            logger.info("Calling image {0}".format(image))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
                working_dir=working_dir,
                parameters=[
                    "perl",
                    "/home/biodocker/bin/NOVOPlasty.pl",
                    "-c",
                    "config.txt",
                ],
            )

            # Write the log, and contigs FASTA files from the local temporary
            # directory into the file store
            log_file_id = utilities.writeGlobalFile(fileStore, log_file_name)
            contigs_file_id = utilities.writeGlobalFile(fileStore, contigs_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            contigs_file_id = None
            log_file_id = None

        # Return file ids and names for export
        novoplasty_rv = {
            "novoplasty_rv": {
                "log_file": {
                    "id": log_file_id,
                    "name": log_file_name,
                },
                "contigs_file": {
                    "id": contigs_file_id,
                    "name": contigs_file_name,
                },
            }
        }
        novoplasty_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(novoplasty_rv))
        return novoplasty_rv


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
        "-w", "--well-spec", default="A01", help="the well specification"
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
        "--seed-path",
        default=os.sep + os.path.join(*cmps),
        help="path to a .fastq file to be used as seed by assembler, where applicable (defaults to using first read)",
    )
    parser.add_argument(
        "--seed-file",
        default=None,
        help="path to a .fastq file to be used as seed by assembler, where applicable (defaults to using first read)",
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

            # Import local seed file into the file store, if needed
            if options.seed_file:
                seed_file_id = utilities.importFile(
                    toil, os.path.join(options.seed_path, options.seed_file)
                )
            else:
                seed_file_id = None

            # Construct and start the NOVOPlasty job
            novoplasty_job = NovoplastyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                options.config_file,
                seed_file_id=seed_file_id,
            )
            novoplasty_rv = toil.start(novoplasty_job)

        else:

            # Restart the NOVOPlasty job
            novoplasty_rv = toil.restart(novoplasty_job)

        # Export all NOVOPlasty output files from the file store
        utilities.exportFiles(
            toil, options.output_directory, novoplasty_rv["novoplasty_rv"]
        )
