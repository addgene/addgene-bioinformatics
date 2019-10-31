from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from SpadesJob import SpadesJob

import utilities


class WellAssemblyJob(Job):
    """
    Runs a SpadesJob followed by an ApcJob.
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 coverage_cutoff, output_directory, *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of file in job store containing FASTQ Illumina short
            left paired reads
        read_two_file_id : toil.fileStore.FileID
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
    Assemble reads and circularize contigs corresponding to a single
    well.

    """
    # Parse FASTQ data path, plate and well specification, coverage
    # cutoff, and output directory, making the output directory if
    # needed
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
                toil, options.data_path, options.plate_spec, [options.well_spec],
                options.source_scheme)

            # Construct and start the well assembly job
            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.coverage_cutoff,
                options.output_directory,
                )
            well_assembly_rv = toil.start(well_assembly_job)

        else:

            # Restart the well assembly job
            well_assembly_rv = toil.restart(well_assembly_job)

        # Export the SPAdes warnings and log files, and contigs FASTA
        # file from the file store
        utilities.exportFiles(
            toil, options.output_directory, well_assembly_rv['spades_rv'])

        # Export the apc output file, and sequence FASTA file from the
        # file store
        utilities.exportFiles(
            toil, options.output_directory, well_assembly_rv['apc_rv'])
