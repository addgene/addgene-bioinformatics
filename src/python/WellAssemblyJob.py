from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from NovoplastyJob import NovoplastyJob
from SpadesJob import SpadesJob
from ShovillJob import ShovillJob

import utilities


class WellAssemblyJob(Job):
    """
    Runs a SpadesJob or ShovillJob followed by an ApcJob, or a
    NovoplastyJob
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 assembler, coverage_cutoff, output_directory,
                 *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of file in file store containing FASTQ Illumina short
            left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of file in file store containing FASTQ Illumina short
            right paired reads
        assembler : str
            name of assembler to run, from ['spades', 'shovill', 'novoplasty']
        coverage_cutoff : str
            read coverage cutoff value
        output_directory : str
            name of directory for output
        """
        super(WellAssemblyJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        if assembler not in ['spades', 'shovill', 'novoplasty']:
            raise Exception("Unexpected assembler")
        self.assembler = assembler
        self.coverage_cutoff = coverage_cutoff
        self.output_directory = output_directory

    def run(self, fileStore):
        """
        Returns
        -------
        dict
            the return values of the SpadesJob and ApcJob
        """
        if self.assembler == 'spades':
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

        elif self.assembler == 'shovill':
            shovill_job = ShovillJob(
                self.read_one_file_id,
                self.read_two_file_id,
                self.output_directory,
            )
            apc_job = ApcJob(
                shovill_job.rv('shovill_rv', 'contigs_file', 'id'),
                parent_rv=shovill_job.rv()
            )
            final_job = self.addChild(
                shovill_job).addChild(
                    apc_job)

        else:  # self.assembler == 'novoplasty'
            novoplasty_job = NovoplastyJob(
                self.read_one_file_id,
                self.read_two_file_id,
            )
            final_job = self.addChild(
                novoplasty_job)

        # Return file ids, and names for export
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
    parser.add_argument('-a', '--assembler', default="spades",
                        choices=['spades', 'shovill', 'novoplaty'],
                        help="name of assembler to run")
    parser.add_argument('-c', '--coverage-cutoff', default="100",
                        help="the coverage cutoff")
    parser.add_argument('-o', '--output-directory', default=None,
                        help="the directory containing all output files")

    # Define and make the output directory, if needed
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

            # Construct and start the well assembly job
            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.assembler,
                options.coverage_cutoff,
                options.output_directory,
                )
            well_assembly_rv = toil.start(well_assembly_job)

        else:

            # Restart the well assembly job
            well_assembly_rv = toil.restart(well_assembly_job)

        # Export the assembler output files from the file store
        utilities.exportFiles(toil, options.output_directory,
                              well_assembly_rv[options.assembler + "_rv"])

        if self.assembler in ['spades', 'shovill']:
            # Export the apc output file, and sequence FASTA file from
            # the file store
            utilities.exportFiles(
                toil, options.output_directory, well_assembly_rv['apc_rv'])
