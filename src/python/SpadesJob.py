import os
import subprocess

from toil.job import Job

import utilities


class SpadesJob(Job):
    """
    Accepts paired-end Illumina reads for assembly with a coverage
    cutoff specification
    """

    def __init__(self, read_one_file_id, read_two_file_id,
                 coverage_cutoff, output_directory, parent_rv={},
                 *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : str
            id of the file in the job store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : str
            id of the file in the job store containing FASTQ Illumina
            short right paired reads
        coverage_cutoff : str
            read coverage cutoff value
        output_directory : str
            name of directory for output
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
        dict
            file ids of warnings and spades log files, and contigs
            FASTA file
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
        warnings_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, "warnings.log")
        spades_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, "spades.log")
        contigs_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, "contigs.fasta")

        # Return file ids for export
        spades_rv = {
            'spades_rv': {
                'warnings_file_id': warnings_file_id,
                'spades_file_id': spades_file_id,
                'contigs_file_id': contigs_file_id,
                }
            }
        spades_rv.update(self.parent_rv)
        return spades_rv
