from argparse import ArgumentParser
import logging
import os
from pathlib import Path

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class SpadesJob(Job):
    """
    Accepts paired-end Illumina reads for assembly with a coverage
    cutoff specification.
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 output_directory, parent_rv={},
                 read_one_file_name="R1.fastq.gz",
                 read_two_file_name="R2.fastq.gz",
                 config_file_path=None,
                 *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short right paired reads
        output_directory : str
            name of directory for output
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(SpadesJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.output_directory = output_directory
        self.parent_rv = parent_rv
        self.config_file_path = config_file_path

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of warnings and spades log files, and
            contigs FASTA file
        """
        # Expected output file names
        warnings_file_name = "warnings.log"
        spades_file_name = "spades.log"
        contigs_file_name = "contigs.fasta"

        try:
            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, self.read_one_file_name)
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, self.read_two_file_name)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call spades.py
            # TODO: Specify the container on construction
            image = "ralatsdio/spades:v3.13.1"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))
            parameters = ["spades.py",
                          "-1",
                          read_one_file_path,
                          "-2",
                          read_two_file_path,
                          "-o",
                          os.path.join(working_dir, self.output_directory),
                         ]

            if self.config_file_path is not None:
                parsed_params = utilities.parseConfigFile(self.config_file_path)
                parameters.extend(parsed_params)
                logger.info("Adding parsed params to CLI call: " + str(parsed_params))

            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
                working_dir=working_dir,
                parameters=parameters)

            # Write the warnings and spades log files, and contigs FASTA
            # file from the local temporary directory into the file store
            warnings_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, warnings_file_name)
            spades_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, spades_file_name)
            contigs_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, contigs_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            warnings_file_id = None
            spades_file_id = None
            contigs_file_id = None

        # Return file ids and names for export
        spades_rv = {
            'spades_rv': {
                'warnings_file': {
                    'id': warnings_file_id,
                    'name': warnings_file_name,
                },
                'spades_file': {
                    'id': spades_file_id,
                    'name': spades_file_name,
                },
                'contigs_file': {
                    'id': contigs_file_id,
                    'name': contigs_file_name,
                },
            }
        }
        spades_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(spades_rv))
        return spades_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification, coverage
    # cutoff, and output directory, making the output directory if
    # needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
    cmps.extend(["dat", "miscellaneous"])
    parser.add_argument('-d', '--data-path',
                        default=os.sep + os.path.join(*cmps),
                        help="path containing plate and well FASTQ source")
    parser.add_argument('-s', '--source-scheme',
                        default="file",
                        help="scheme used for the source URL")
    parser.add_argument('-p', '--plate-spec', default="A11967A_sW0154",
                        help="the plate specification")
    parser.add_argument('-w', '--well-spec', default="B01",
                        help="the well specification")
    parser.add_argument('-o', '--output-directory', default=None,
                        help="the directory containing all output files")
    parser.add_argument("-c", "--config", default=None,
                        help="a .ini file with args to be passed to SPAdes")  

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
                toil, options.data_path, options.plate_spec,
                [options.well_spec], options.source_scheme)

            # Construct and start the SPAdes job
            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.output_directory,
                config_file_path=str(Path(options.config).absolute()) if options.config is not None else None,
                )
            spades_rv = toil.start(spades_job)

        else:

            # Restart the SPAdes job
            spades_rv = toil.restart(spades_job)

        # Export all SPAdes output files from the file store
        utilities.exportFiles(
            toil, options.output_directory, spades_rv['spades_rv']
        )