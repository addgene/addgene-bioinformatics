import subprocess

from toil.job import Job

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
        dict
            file ids of output and alignment files, and circular
            sequence FASTA file
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

        # Return file ids for export
        apc_rv = {
            'apc_rv': {
                'output_file_id': output_file_id,
                'sequence_file_id': sequence_file_id,
                }
            }
        apc_rv.update(self.parent_rv)
        return apc_rv
