from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities


class MasurcaJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using MaSuRCA.
    """
    def __init__(self, read_one_file_id, read_two_file_id,
                 parent_rv={},
                 read_one_file_name="R1.fastq.gz",
                 read_two_file_name="R2.fastq.gz",
                 threads="1",
                 *args, **kwargs):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short right paired reads
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(MasurcaJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.threads = threads
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of log and corrections files, and
            contigs FASTA file
        """
        # Expected output file names
        ca0_file_name = "runCA0.out"  # CABOG stdout for initial stages: store building, overlapping, unitigging
        ca1_file_name = "runCA1.out"  # CABOG stdout for unitig consensus
        ca2_file_name = "runCA2.out"  # CABOG stdout for scaffolder
        sr1_file_name = "super1.err"  # STDERR for super reads code that generates super reads from PE reads
        ecl_file_name = "error_correct.log"  # Log file for error correction
        scaffolds_file_name = "final.genome.scf.fasta"

        try:
            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, self.read_one_file_name)
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, self.read_two_file_name)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call MaSuRCA
            # TODO: Specify the container on construction
            working_dir = fileStore.localTempDir
            apiDockerCall(
                self,
                image='ralatsdio/masurca:v1.0.9',
                volumes={working_dir: {'bind': working_dir, 'mode': 'rw'}},
                working_dir=working_dir,
                parameters=["assemble.sh",
                            "--short1",
                            read_one_file_path,
                            "--short2",
                            read_two_file_path,
                            "--threads",
                            self.threads,
                            ])

            # Write the MaSuRCA log and corrections files, and contigs
            # FASTA file from the local temporary directory into the file
            # store
            ca0_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, ca0_file_name)
            ca1_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, ca1_file_name)
            ca2_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, ca2_file_name)
            sr1_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, sr1_file_name)
            ecl_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, ecl_file_name)
            scaffolds_file_id = utilities.writeGlobalFile(
                fileStore, self.output_directory, scaffolds_file_name)

        except Exception as exc:
            # Ensure expectred return values on exceptions
            ca0_file_id = None
            ca1_file_id = None
            ca2_file_id = None
            sr1_file_id = None
            ecl_file_id = None
            scaffolds_file_id = None

        # Return file ids and names for export
        masurca_rv = {
            'masurca_rv': {
                'ca0_file': {
                    'id': ca0_file_id,
                    'name': ca0_file_name,
                },
                'ca1_file': {
                    'id': ca1_file_id,
                    'name': ca1_file_name,
                },
                'ca2_file': {
                    'id': ca2_file_id,
                    'name': ca2_file_name,
                },
                'sr1_file': {
                    'id': sr1_file_id,
                    'name': sr1_file_name,
                },
                'ecl_file': {
                    'id': ecl_file_id,
                    'name': ecl_file_name,
                },
                'scaffolds_file': {
                    'id': scaffolds_file_id,
                    'name': scaffolds_file_name,
                },
            }
        }
        masurca_rv.update(self.parent_rv)
        return masurca_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification, and output
    # directory, making the output directory if needed
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

            # STARTHERE: Sort how to handle assemble.cfg or assemble.sh

            # Import the local read files into the file store
            read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
                toil, options.data_path, options.plate_spec,
                [options.well_spec], options.source_scheme)

            # Construct and start the MaSuRCA job
            masurca_job = MasurcaJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.output_directory,
                )
            masurca_rv = toil.start(masurca_job)

        else:

            # Restart the MaSuRCA job
            masurca_rv = toil.restart(masurca_job)

        # Export the MaSuRCA log and corrections files, and contigs
        # FASTA file from the file store
        utilities.exportFiles(
            toil, options.output_directory, masurca_rv['masurca_rv'])
