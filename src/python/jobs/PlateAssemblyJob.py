from argparse import ArgumentParser
from glob import glob
import logging
import os
from pathlib import Path
import re

import boto
from toil.job import Job
from toil.common import Toil

from WellAssemblyJob import WellAssemblyJob

import utilities

logger = logging.getLogger(__name__)


class PlateAssemblyJob(Job):
    """
    Runs a WellAssemblyJob for each well of a plate for which both
    FASTQ read files were found.
    """

    def __init__(
        self,
        well_specs,
        read_one_file_ids,
        read_two_file_ids,
        plate_spec,
        assembler,
        config_file_id,
        config_file_name,
        *args,
        **kwargs
    ):
        """
        Parameters
        ----------
        well_specs : list of str
            specification for each well for which both FASTQ read
            files were found
        read_one_file_ids : list of toil.fileStore.FileID
            ids of files in file store containing FASTQ Illumina short
            left paired reads
        read_two_file_ids : list of toil.fileStore.FileID
            ids of files in file store containing FASTQ Illumina short
            right paired reads
        plate_spec : str
            specification for plate containing the specified wells
        assembler : str
            name of assembler to run, from utilities.ASSEMBLERS_TO_RUN
        config_file_id : toil.fileStore.FileID
            id of the file in the file store containing assembler args
        config_file_name : str
            name of the file in the file store containing assembler args
        """
        super(PlateAssemblyJob, self).__init__(*args, **kwargs)
        self.well_specs = well_specs
        self.read_one_file_ids = read_one_file_ids
        self.read_two_file_ids = read_two_file_ids
        self.plate_spec = plate_spec
        if assembler not in utilities.ASSEMBLERS_TO_RUN:
            raise Exception("Unexpected assembler")
        self.assembler = assembler
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name

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
            logger.info(
                "Creating well assembly job {0}".format(
                    self.plate_spec + "_" + self.well_specs[iW]
                )
            )
            # Note that exceptions are caught in the well assembly job
            well_assembly_rvs.append(
                self.addChild(
                    WellAssemblyJob(
                        self.read_one_file_ids[iW],
                        self.read_two_file_ids[iW],
                        self.assembler,
                        self.config_file_id,
                        self.config_file_name,
                        self.plate_spec + "_" + self.well_specs[iW],
                    )
                ).rv()
            )
        return well_assembly_rvs


if __name__ == "__main__":
    """
    Assemble reads for each well of a specified plate for which both
    FASTQ read files are found.
    """
    # Parse FASTQ data path, plate specification, coverage cutoff, and
    # output directory
    parser = ArgumentParser()
    Job.Runner.addToilOptions(parser)
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
    cmps.extend(["dat", "miscellaneous"])
    parser.add_argument(
        "-d",
        "--data-path",
        default=os.sep + os.path.join(*cmps),
        help="path containing plate and well FASTQ source",
    )
    parser.add_argument(
        "-s", "--source-scheme", default="file", help="scheme used for the source URL"
    )
    parser.add_argument(
        "-l", "--plate-spec", default="A11967B_sW0154", help="the plate specification"
    )
    parser.add_argument(
        "-a",
        "--assembler",
        default="spades",
        choices=utilities.ASSEMBLERS_TO_RUN,
        help="name of the assembler to run",
    )
    cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-1]
    parser.add_argument(
        "-c",
        "--config-path",
        default=os.sep + os.path.join(*cmps),
        help="path to a .ini file with args to be passed to the assembler",
    )
    parser.add_argument(
        "-f",
        "--config-file",
        default="Assembler.ini",
        help="path to a .ini file with args to be passed to the assembler",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        default=None,
        help="the directory containing all output files",
    )

    # Define and make the output directory, if needed
    options = parser.parse_args()
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
            if options.source_scheme == "file":

                # Find read one and two files
                read_one_files = glob(
                    os.path.join(
                        options.data_path,
                        "{0}_FASTQ".format(options.plate_spec),
                        "{0}_*_R1_001.fastq.gz".format(options.plate_spec),
                    )
                )
                read_two_files = glob(
                    os.path.join(
                        options.data_path,
                        "{0}_FASTQ".format(options.plate_spec),
                        "{0}_*_R2_001.fastq.gz".format(options.plate_spec),
                    )
                )

            elif options.source_scheme == "s3":

                # Find the bucket name and path to the data source
                cmps = options.data_path.split(os.sep)
                bucket_name = cmps[0]
                bucket_cmps = cmps[1:]
                bucket_cmps.append("{0}_FASTQ".format(options.plate_spec))
                bucket_path = os.path.join(*bucket_cmps)

                # Find read one and two files
                s3 = boto.connect_s3()
                bucket = s3.get_bucket(bucket_name)
                read_one_files = [
                    f.name for f in bucket.list(bucket_path) if re.search("R1", f.name)
                ]
                read_two_files = [
                    f.name for f in bucket.list(bucket_path) if re.search("R2", f.name)
                ]

            # Ensure a read two file exists for each read one file
            for read_one_file in read_one_files:
                read_two_file = read_one_file.replace("R1", "R2")
                if read_two_file in read_two_files:
                    well_specs.append(p.search(read_one_file).group(1))

            # Import local read file pairs
            read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
                toil,
                options.data_path,
                options.plate_spec,
                well_specs,
                options.source_scheme,
            )

            # Import local config file into the file store
            config_file_id = utilities.importConfigFile(
                toil, os.path.join(options.config_path, options.config_file)
            )

            # Construct and start the plate assembly job
            plate_assembly_job = PlateAssemblyJob(
                well_specs,
                read_one_file_ids,
                read_two_file_ids,
                options.plate_spec,
                options.assembler,
                config_file_id,
                options.config_file,
            )
            well_assembly_rvs = toil.start(plate_assembly_job)

        else:

            # Restart the well assembly job
            well_assembly_rvs = toil.restart(plate_assembly_job)

        # Export all assembler output files, and apc output files, if
        # apc is required by the asembler, from the file store
        utilities.exportWellAssemblyFiles(
            toil,
            options.assembler,
            options.output_directory,
            well_specs,
            well_assembly_rvs,
        )
