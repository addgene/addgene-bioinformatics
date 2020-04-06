from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class SkesaJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using SKESA.
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 parent_rv={},
                 read_one_file_name="R1.fastq.gz",
                 read_two_file_name="R2.fastq.gz",
                 output_file="contigs.fa",
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
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(SkesaJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.output_file = output_file
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of log and corrections files, and
            contigs FASTA file
        """
        # Expected output file names
        # TODO: Find a way to capture stderr
        # log_file_name = "skesa.log"
        contigs_file_name = self.output_file

        try:
            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, self.read_one_file_name)
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, self.read_two_file_name)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call SKESA
            # TODO: Specify the container on construction
            image = "ralatsdio/skesa:v2.3.0"
            working_dir = fileStore.localTempDir
            logger.info("Calling image {0}".format(image))
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
                working_dir=working_dir,
                parameters=["skesa",
                            "--fastq",
                            ",".join([read_one_file_path, read_two_file_path]),
                            "--contigs_out",
                            os.path.join(working_dir, contigs_file_name),
                            ],
                # stderr=True,
                # steamfile=log_file_name,
                )

            # Write the SKESA log, and contigs FASTA file from the
            # local temporary directory into the file store
            # log_file_id = utilities.writeGlobalFile(
            #     fileStore, log_file_name)
            contigs_file_id = utilities.writeGlobalFile(
                fileStore, contigs_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            # log_file_id = None
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            contigs_file_id = None

        # Return file ids and names for export
        skesa_rv = {
            'skesa_rv': {
                # 'log_file': {
                #     'id': log_file_id,
                #     'name': log_file_name,
                # },
                'contigs_file': {
                    'id': contigs_file_id,
                    'name': contigs_file_name,
                },
            }
        }
        skesa_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(skesa_rv))
        return skesa_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification, and output
    # directory, making the output directory if needed
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

            # Construct and start the SKESA job
            skesa_job = SkesaJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                )
            skesa_rv = toil.start(skesa_job)

        else:

            # Restart the SKESA job
            skesa_rv = toil.restart(skesa_job)

        # Export all SKESA output files from the file store
        utilities.exportFiles(
            toil, options.output_directory, skesa_rv['skesa_rv']
        )
