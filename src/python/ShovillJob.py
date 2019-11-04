from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities


class ShovillJob(Job):
    """
    Accepts paired-end Illumina reads for assembly with SPAdes, SKESA,
    Megahit, or Velvet performing desirable pre- and post-processing.
    """

    def __init__(self, read_one_file_id, read_two_file_id,
                 output_directory, parent_rv={},
                 read_one_file_name="R1.fastq.gz", read_two_file_name="R2.fastq.gz",
                 assembler="spades",  # Use spades, skesa, megahit, or velvet
                 cpus=16,             # Use this number of CPUs (0=ALL)
                 ram=4.00,            # Try to keep RAM usage below this many GB
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
        output_directory : str
            name of directory for output
        parent_rv : dict
            dictionary of return values from the parent job
        assembler : str
            One of skesa, velvet, megahit, or spades
        cpus : int
            Number of CPUs to use (0=ALL)
        ram : flt
            RAM usage limit [GB]
        """
        super(ShovillJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.output_directory = output_directory
        self.parent_rv = parent_rv
        self.assembler = assembler
        self.cpus = cpus
        self.ram = ram

    def run(self, fileStore):
        """
        Returns
        -------
        TODO: Complete when it is clear what output files are created
        dict of toil.fileStore.FileID and str
            file ids and names of assorted output files.
        """
        # Read the read files from the file store into the local
        # temporary directory
        read_one_file_path = utilities.readGlobalFile(
            fileStore, self.read_one_file_id, self.read_one_file_name)
        read_two_file_path = utilities.readGlobalFile(
            fileStore, self.read_two_file_id, self.read_two_file_name)

        # Mount the Toil local temporary directory to the same path in
        # the container, and use the path as the working directory in
        # the container, then call shovill
        # TODO: Specify the container on construction
        working_dir = fileStore.localTempDir
        apiDockerCall(
            self,
            image='bioboxes/shovill',
            volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
            working_dir=working_dir,
            parameters=["shovill",
                        "-1",
                        read_one_file_path,
                        "-2",
                        read_two_file_path,
                        "--outdir",
                        os.path.join(working_dir, self.output_directory),
                        "--cpus",
                        self.cpus,
                        "--ram",
                        self.ram,
                        "--assembler",
                        self.assembler,
                        ])

        # TODO: >>> Start here
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
    cmps = str(os.path.abspath(__file__)).split(os.sep)
    parser.add_argument('-d', '--data-directory',
                        default=os.sep + os.path.join(*cmps[0:-3], "dat", "miscellaneous"),
                        help="the directory containing FASTQ read data files")
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
                toil, options.data_directory, options.plate_spec, [options.well_spec])

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
