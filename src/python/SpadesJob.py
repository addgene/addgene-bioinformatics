from argparse import ArgumentParser
import os
import subprocess

from toil.job import Job
from toil.common import Toil

import utilities


class SpadesJob(Job):
    """
    Accepts paired-end Illumina reads for assembly with a coverage
    cutoff specification.
    """

    def __init__(self, read_one_file_id, read_two_file_id,
                 coverage_cutoff, output_directory, parent_rv={},
                 *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of the file in the job store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of the file in the job store containing FASTQ Illumina
            short right paired reads
        coverage_cutoff : str
            read coverage cutoff value
        output_directory : str
            name of directory for output
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(SpadesJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.coverage_cutoff = coverage_cutoff
        self.output_directory = output_directory
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of warnings and spades log files, and
            contigs FASTA file
        """
        # Read the read files from the file store into the local
        # temporary directory
        read_one_file_path = utilities.readGlobalFile(
            fileStore, self.read_one_file_id, "R1.fastq.gz")
        read_two_file_path = utilities.readGlobalFile(
            fileStore, self.read_two_file_id, "R2.fastq.gz")

        # Call spades.py
        subprocess.check_call([
            "spades.py",
            "-1",
            read_one_file_path,
            "-2",
            read_two_file_path,
            "-o",
            os.path.join(fileStore.localTempDir, self.output_directory),
            "--cov-cutoff",
            self.coverage_cutoff,
            ])

        # Write the warnings and spades log files, and contigs FASTA
        # file from the local temporary directory into the file store
        warnings_file_name = "warnings.log"
        warnings_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, warnings_file_name)
        spades_file_name = "spades.log"
        spades_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, spades_file_name)
        contigs_file_name = "contigs.fasta"
        contigs_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, contigs_file_name)

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
        return spades_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.

    """
    # Parse FASTQ data directory, plate and well specification,
    # coverage cutoff, and output directory, making the output
    # directory if needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('-d', '--data-directory',
                        default=os.path.join("..", "..", "dat", "miscellaneous"),
                        help="the directory containing FASTQ read data files")
    parser.add_argument('-p', '--plate-spec', default="A11967A_sW0154",
                        help="the plate specification")
    parser.add_argument('-w', '--well-spec', default="A01",
                        help="the well specification")
    parser.add_argument('-c', '--coverage-cutoff', default="100",
                        help="the coverage cutoff")
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
                toil, options.data_directory,
                options.plate_spec, [options.well_spec])

            # Construct and start the SPAdes job
            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.coverage_cutoff,
                options.output_directory,
                )
            spades_rv = toil.start(spades_job)

        else:

            # Restart the SPAdes job
            spades_rv = toil.restart(spades_job)

        # Export the SPAdes warnings and log files, and contigs FASTA
        # file from the file store
        utilities.exportFiles(
            toil, options.output_directory, spades_rv['spades_rv'])
