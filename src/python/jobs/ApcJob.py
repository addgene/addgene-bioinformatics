from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class ApcJob(Job):
    """
    Tests each sequence in a supplied fasta file for self-overlap. If
    overlap is found, the 5' copy of the overlapping sequence is
    trimmed, the ends joined, and the join moved to the middle of the
    output sequence. The join location is appended to the sequence
    header, to allow confirmation.
    """

    def __init__(
        self,
        contigs_file_id,
        parent_rv={},
        contigs_file_name="contigs.fasta",
        base_file_name="apc",
        *args,
        **kwargs
    ):
        """
        Parameters
        ----------
        contigs_file_id : str
            id of the file in file store containing FASTA contigs
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(ApcJob, self).__init__(*args, **kwargs)
        self.contigs_file_id = contigs_file_id
        self.contigs_file_name = contigs_file_name
        self.base_file_name = "apc"
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of output and alignment files, and
            circular sequence FASTA file
        """
        # Expected output file names
        log_file_name = self.base_file_name + ".log"
        sequence_file_name = self.base_file_name + ".1.fa"

        try:
            # Read the contigs FASTA file from the file store into the
            # local temporary directory
            utilities.readGlobalFile(
                fileStore, self.contigs_file_id, self.contigs_file_name
            )

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call apc.pl
            # TODO: Specify the container on construction
            image = "ralatsdio/apc:v0.1.0"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
                working_dir=working_dir,
                parameters=[
                    "apc.pl",
                    "-b",
                    self.base_file_name,
                    self.contigs_file_name,
                ],
            )

            # Write the output file from the local temporary directory
            # into the file store
            log_file_id = utilities.writeGlobalFile(fileStore, log_file_name)
            sequence_file_id = utilities.writeGlobalFile(fileStore, sequence_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            log_file_id = None
            sequence_file_id = None

        # Return file ids and names for export
        apc_rv = {
            "apc_rv": {
                "log_file": {"id": log_file_id, "name": log_file_name,},
                "sequence_file": {"id": sequence_file_id, "name": sequence_file_name,},
            }
        }
        apc_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(apc_rv))
        return apc_rv


if __name__ == "__main__":
    """
    Circularize contigs corresponding to a single well.
    """
    # Parse FASTA data path and file name, and output directory,
    # making the output directory if needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-1]
    cmps.extend(["A11967A_sW0154_B01"])
    parser.add_argument(
        "-d",
        "--data-path",
        default=os.sep + os.path.join(*cmps),
        help="path containing the FASTA source",
    )
    parser.add_argument(
        "-s", "--source-scheme", default="file", help="scheme used for the source URL"
    )
    parser.add_argument(
        "-f", "--file-name", default="contigs.fasta", help="the FASTA contigs file"
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        default=None,
        help="the directory containing all output files",
    )
    options = parser.parse_args()
    if options.output_directory is None:
        options.output_directory = os.path.relpath(options.data_path)
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)

    # Work within the Toil context manager
    with Toil(options) as toil:
        if not toil.options.restart:

            # Import the local contigs file into the file store
            contigs_file_id = utilities.importFile(
                toil,
                os.path.join(options.data_path, options.file_name),
                options.source_scheme,
            )

            # Construct and start the APC job
            apc_job = ApcJob(contigs_file_id,)
            apc_rv = toil.start(apc_job)

        else:

            # Restart the APC job
            apc_rv = toil.restart(apc_job)

        # Export all apc output files from the file store
        utilities.exportFiles(toil, options.output_directory, apc_rv["apc_rv"])
