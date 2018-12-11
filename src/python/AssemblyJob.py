from toil.job import Job

from ApcJob import ApcJob
from SpadesJob import SpadesJob


class AssemblyJob(Job):
    """
    Runs a SpadesJob followed by an ApcJob
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 coverage_cutoff, output_directory, *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : str
            id of file in job store containing FASTQ Illumina short
            left paired reads
        read_two_file_id : str
            id of file in job store containing FASTQ Illumina short
            right paired reads
        coverage_cutoff : str
            read coverage cutoff value
        output_directory : str
            name of directory for output
        """
        super(AssemblyJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.coverage_cutoff = coverage_cutoff
        self.output_directory = output_directory

    def run(self, fileStore):
        """
        Returns
        -------
        dict
            the return values of the SpadesJob and ApcJob
        """
        spades_job = SpadesJob(
            self.read_one_file_id,
            self.read_two_file_id,
            self.coverage_cutoff,
            self.output_directory,
            )
        apc_job = ApcJob(
            spades_job.rv('spades_rv', 'contigs_file_id'),
            parent_rv=spades_job.rv()
            )
        final_job = self.addChild(
            spades_job) # .addChild(
        # apc_job)

        # Return file ids for export
        return final_job.rv()
