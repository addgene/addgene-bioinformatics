from argparse import ArgumentParser
from glob import glob
import os
import re

from toil.job import Job
from toil.common import Toil

from WellAssemblyJob import WellAssemblyJob

import utilities


class PlateAssemblyJob(Job):
    """
    Runs a WellAssemblyJob for each well of a plate for which both
    FASTQ read files were found
    """
    def __init__(self, well_specs, read_one_file_ids, read_two_file_ids,
                 plate_spec, coverage_cutoff, *args, **kwargs):
        """
        Parameters
        ----------
        well_specs : list of str
            Specification for each well for which both FASTQ read
            files were found
        read_one_file_ids : list of toil.fileStore.FileID
            ids of files in job store containing FASTQ Illumina short
            left paired reads
        read_two_file_ids : list of toil.fileStore.FileID
            ids of files in job store containing FASTQ Illumina short
            right paired reads
        plate_spec : str
            Specification for plate containing the specified wells
        coverage_cutoff : str
            read coverage cutoff value
        """
        super(PlateAssemblyJob, self).__init__(*args, **kwargs)
        self.well_specs = well_specs
        self.read_one_file_ids = read_one_file_ids
        self.read_two_file_ids = read_two_file_ids
        self.plate_spec = plate_spec
        self.coverage_cutoff = coverage_cutoff

    def run(self, fileStore):
        """
        Returns
        -------
        dict
            the return values of the WellAssemblyJobs
        """
        well_assembly_rvs = []
        nW = len(self.well_specs)
        for iW in range(nW):
            well_assembly_rvs.append(
                self.addChild(
                    WellAssemblyJob(
                        self.read_one_file_ids[iW],
                        self.read_two_file_ids[iW],
                        self.coverage_cutoff,
                        self.plate_spec + "_" + self.well_specs[iW],
                        cores=1,
                        disk="2G",
                        memory="3G",
                        )).rv())
        return well_assembly_rvs


if __name__ == "__main__":
    """
    Assemble reads for each well of a specified plate for which both
    FASTQ read files are found
    """

    # Parse FASTQ data directory, plate specification, coverage
    # cutoff, and output directory
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    parser.add_argument('-d', '--data-directory',
                        default="../../dat/miscellaneous",
                        help="the directory containing FASTQ read data files")
    parser.add_argument('-p', '--plate-spec', default="A11967B_sW0154",
                        help="the plate specification")
    parser.add_argument('-c', '--coverage-cutoff', default="100",
                        help="the coverage cutoff")
    parser.add_argument('-o', '--output-directory', default=None,
                        help="the directory containing all output files")
    options = parser.parse_args()

    # Define and make the output directory, if needed
    if options.output_directory is None:
        options.output_directory = options.plate_spec
    if not os.path.exists(options.output_directory):
        os.mkdir(options.output_directory)

    # Work within the Toil context manager
    with Toil(options) as toil:
        if not toil.options.restart:

            # Import local read file pairs into the file store
            p = re.compile(r"{0}_(\w+)_R1_001".format(options.plate_spec))
            well_specs = []
            read_one_file_ids = []
            read_two_file_ids = []
            read_one_files = glob(os.path.join(
                options.data_directory,
                "{0}_FASTQ/{0}_*_R1_001.fastq.gz".format(options.plate_spec)))
            for read_one_file in read_one_files:
                read_two_file = read_one_file.replace("R1", "R2")
                if os.path.exists(read_two_file):
                    well_specs.append(p.search(read_one_file).group(1))
            read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
                toil, options.data_directory, options.plate_spec, well_specs)

            # Construct and start the plate assembly job
            plate_assembly_job = PlateAssemblyJob(
                well_specs,
                read_one_file_ids,
                read_two_file_ids,
                options.plate_spec,
                options.coverage_cutoff,
                cores=2,
                disk="3G",
                memory="4G",
                )
            well_assembly_rvs = toil.start(plate_assembly_job)

        else:

            # Restart the well assembly job
            well_assembly_rvs = toil.restart(plate_assembly_job)

        # Export needed files created by each well assembly job
        nW = len(well_specs)
        for iW in range(nW):

            # Define and make the well output directory, if needed
            well_output_directory = os.path.join(
                options.output_directory, well_specs[iW])
            if not os.path.exists(well_output_directory):
                os.mkdir(well_output_directory)

            # Export the SPAdes warnings and log files, and contigs
            # FASTA file from the file store
            utilities.exportFiles(
                toil, well_output_directory, well_assembly_rvs[iW]['spades_rv'])

            # Export the apc output file, and sequence FASTA file from
            # the file store
            utilities.exportFiles(
                toil, well_output_directory, well_assembly_rvs[iW]['apc_rv'])
