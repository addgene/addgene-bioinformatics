from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class UnicyclerJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using Unicycler.
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 output_directory, parent_rv={},
                 read_one_file_name="R1.fastq.gz",
                 read_two_file_name="R2.fastq.gz",
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
        super(UnicyclerJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.output_directory = output_directory
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of log, and contigs FASTA file
        """
        # Expected output file names
        log_file_name = "unicycler.log"
        contigs_file_name = "assembly.fasta"
        graph_file_name = "assembly.gfa"

        try:
            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, self.read_one_file_name)
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, self.read_two_file_name)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call Unicycler
            # TODO: Specify the container on construction
            image = "ralatsdio/unicycler:v0.4.7"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
                working_dir=working_dir,
                parameters=["unicycler",
                            "--short1",
                            read_one_file_path,
                            "--short2",
                            read_two_file_path,
                            "--out",
                            os.path.join(working_dir, self.output_directory),
                            ])

            # Write the Unicycler log, and final assembly contigs
            # FASTA and graph GFA v1 files from the local temporary
            # directory into the file store
            log_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, log_file_name)
            contigs_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, contigs_file_name)
            graph_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, graph_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            log_file_id = None
            contigs_file_id = None
            graph_file_id = None

        # Return file ids and names for export
        unicycler_rv = {
            'unicycler_rv': {
                'log_file': {
                    'id': log_file_id,
                    'name': log_file_name,
                },
                'contigs_file': {
                    'id': contigs_file_id,
                    'name': contigs_file_name,
                },
                'graph_file': {
                    'id': graph_file_id,
                    'name': graph_file_name,
                },
            }
        }
        unicycler_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(unicycler_rv))
        return unicycler_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification, and output
    # directory, making the output directory if needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-3]
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

            # Construct and start the Unicycker job
            unicycler_job = UnicyclerJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.output_directory,
                )
            unicycler_rv = toil.start(unicycler_job)

        else:

            # Restart the Unicycler job
            unicycler_rv = toil.restart(unicycler_job)
            unicycler_rv = toil.restart(unicycler_job)

        # Export the Unicycler log and corrections files, and contigs
        # FASTA file from the file store
        utilities.exportFiles(
            toil, options.output_directory, unicycler_rv['unicycler_rv'])
