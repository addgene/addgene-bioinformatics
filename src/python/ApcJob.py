from argparse import ArgumentParser
import os
import subprocess

from toil.job import Job
from toil.common import Toil

import utilities


class ApcJob(Job):
    """
    Tests each sequence in a supplied fasta file for self-overlap. If
    overlap is found, the 5' copy of the overlapping sequence is
    trimmed, the ends joined, and the join moved to the middle of the
    output sequence. The join location is appended to the sequence
    header, to allow confirmation.
    """

    def __init__(self, contigs_file_id, parent_rv={}, *args, **kwargs):
        """
        Parameters
        ----------
        contigs_file_id : str
            id of the file in job store containing FASTA contigs
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(ApcJob, self).__init__(*args, **kwargs)
        self.contigs_file_id = contigs_file_id
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of output and alignment files, and
            circular sequence FASTA file
        """
        # Read the contigs FASTA file from the file store into the
        # local temporary directory
        contigs_file_name = "contigs.fasta"
        utilities.readGlobalFile(
            fileStore, self.contigs_file_id, contigs_file_name)

        # Call apc.pl
        base_file_name = "apc"
        output_file_name = base_file_name + ".out"
        sequence_file_name = base_file_name + ".1.fa"
        with open(output_file_name, 'w') as output_file:
            subprocess.check_call([
                "apc.pl",
                "-b",
                base_file_name,
                contigs_file_name,
                ], stdout=output_file)

        # Write the output file from the local temporary directory
        # into the file store
        output_file_id = utilities.writeGlobalFile(
            fileStore, output_file_name)
        sequence_file_id = utilities.writeGlobalFile(
            fileStore, sequence_file_name)

        # Return file ids and names for export
        apc_rv = {
            'apc_rv': {
                'output_file': {
                    'id': output_file_id,
                    'name': output_file_name,
                },
                'sequence_file': {
                    'id': sequence_file_id,
                    'name': sequence_file_name,
                },
            }
        }
        apc_rv.update(self.parent_rv)
        return apc_rv



if __name__ == "__main__":
    """
    Circularize contigs corresponding to a single well.

    """
    # Parse FASTA data directory and file name, and output directory,
    # making the output directory if needed
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('-d', '--data-directory',
                        default=os.path.join("..", "..", "dat", "miscellaneous", "A11967A_sW0154_FASTA"),
                        help="the directory containing FASTA contigs files")
    parser.add_argument('-f', '--file-name', default="A01-A11967A-contigs.fasta",
                        help="the FASTA contigs file")
    parser.add_argument('-o', '--output-directory', default=None,
                        help="the directory containing all output files")
    options = parser.parse_args()
    if options.output_directory is None:
        options.output_directory = options.data_directory
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)

    # Work within the Toil context manager
    with Toil(options) as toil:
        if not toil.options.restart:

            # Import the local contigs file into the file store
            contigs_file_id = utilities.importContigsFile(
                toil, options.data_directory, options.file_name)

            # Construct and start the APC job
            apc_job = ApcJob(
                contigs_file_id,
                )
            apc_rv = toil.start(apc_job)

        else:

            # Restart the APC job
            apc_rv = toil.restart(apc_job)

        # Export the apc output file, and sequence FASTA file from the
        # file store
        utilities.exportFiles(
            toil, options.output_directory, apc_rv['apc_rv'])
