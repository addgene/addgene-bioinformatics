from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from SpadesJob import SpadesJob

import utilities


class WellAssemblyJob(Job):
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
        super(WellAssemblyJob, self).__init__(*args, **kwargs)
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
            spades_job.rv('spades_rv', 'contigs_file', 'id'),
            parent_rv=spades_job.rv()
            )
        final_job = self.addChild(
            spades_job).addChild(
                apc_job)

        # Return file ids for export
        return final_job.rv()


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """

    # Parese data directory, plate and well specification, and
    # coverage cutoff
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('-d', '--data-directory',
                        default="../../dat/miscellaneous",
                        help="the plate specification")
    parser.add_argument('-p', '--plate-spec', default="A11967A_sW0154",
                        help="the plate specification")
    parser.add_argument('-w', '--well-spec', default="A01",
                        help="the well specification")
    parser.add_argument('-c', '--coverage-cutoff', default="100",
                        help="the coverage cutoff")
    options = parser.parse_args()

    # Work within the Toil context manager
    with Toil(options) as toil:
        if not toil.options.restart:

            # Import the local read files into the file store
            read_one_file_id = utilities.importFile(
                toil, os.path.join("../../dat/miscellaneous",
                                   "{0}_FASTQ/{0}_{1}_R1_001.fastq.gz".format(
                                       options.plate_spec, options.well_spec)))
            read_two_file_id = utilities.importFile(
                toil, os.path.join("../../dat/miscellaneous",
                                   "{0}_FASTQ/{0}_{1}_R2_001.fastq.gz".format(
                                       options.plate_spec, options.well_spec)))

            # Construct and start the well assembly job
            assembly_job = WellAssemblyJob(
                read_one_file_id,
                read_two_file_id,
                options.coverage_cutoff,
                options.plate_spec,
                )
            assembly_rv = toil.start(assembly_job)

        else:

            # Restart the well assembly job
            assembly_rv = toil.restart(assembly_job)

        # Export the SPAdes warnings and log files, and contigs FASTA
        # file from the file store
        utilities.exportFiles(toil, assembly_rv['spades_rv'])

        # Export the apc output file, and sequence FASTA file from the
        # file store
        utilities.exportFiles(toil, assembly_rv['apc_rv'])
