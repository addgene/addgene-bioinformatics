from argparse import ArgumentParser
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities


class NovoplastyJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using NOVOPlasty.
    """

    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        output_directory,
        parent_rv={},
        read_one_file_name="R1.fastq.gz",
        read_two_file_name="R2.fastq.gz",
        *args,
        **kwargs
        ):
        """
        Parameters
        ----------
        read_one_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short left paired reads
        read_two_file_id : toil.fileStore.FileID
            id of the file in the file store containing FASTQ Illumina
            short right paired reads
        output_directory : str
            name of directory for output
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(NovoplastyJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_one_file_name = read_one_file_name
        self.read_two_file_id = read_two_file_id
        self.read_two_file_name = read_two_file_name
        self.output_directory = output_directory
        self.parent_rv = parent_rv

    def run(self, fileStore):
        """
        Returns
        -------
        dict of toil.fileStore.FileID and str
            file ids and names of log and corrections files, and
            contigs FASTA file
        """
        # Read the read files from the file store into the local
        # temporary directory
        read_one_file_path = utilities.readGlobalFile(
            fileStore, self.read_one_file_id, self.read_one_file_name
        )
        read_two_file_path = utilities.readGlobalFile(
            fileStore, self.read_two_file_id, self.read_two_file_name
        )

        # TODO: Select a read sequence as the seed here, and write it
        # into the local temporary directory

        # Write the NOVOPlasty config file into the local temporary
        # directory
        working_dir = fileStore.localTempDir
        with open(os.path.join(working_dir, "config.txt"), "w+") as f:
            config = """Project:
-----------------------
Project name          = Assembly
Type                  = mito
Genome Range          = 1800-35000
K-mer                 = 121
Max memory            = 3
Extended log          =
Save assembled reads  =
Seed Input            = Seed.fasta
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = 251
Insert size           = 500
Platform              = illumina
Single/Paired         = PE
Combined reads        = 
Forward reads         = {read_one_file_path}
Reverse reads         = {read_two_file_path}

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3,
Use Quality Scores    = no
            """.format(
                read_one_file_path=read_one_file_path,
                read_two_file_path=read_two_file_path,
            )
            f.write(config)

        # Mount the Toil local temporary directory to the same path in
        # the container, and use the path as the working directory in
        # the container, then call NOVOPlasty
        # TODO: Specify the container on construction
        apiDockerCall(
            self,
            image="ralatsdio/novoplasty:v3.7.0",
            volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
            working_dir=working_dir,
            parameters=["perl",
                        "/home/biodocker/bin/NOVOPlasty.pl",
                        "-c",
                        "config.txt",
                        ],
        )

        # Write the contigs FASTA file from the local temporary
        # directory into the file store
        contigs_file_name = "contigs.fa"
        contigs_file_id = utilities.writeGlobalFile(
            fileStore, self.output_directory, contigs_file_name
        )

        # Return file ids and names for export
        novoplasty_rv = {
            "novoplasty_rv": {
                "contigs_file": {"id": contigs_file_id, "name": contigs_file_name},
            }
        }
        novoplasty_rv.update(self.parent_rv)
        return novoplasty_rv


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
        "-p", "--plate-spec", default="A11967A_sW0154", help="the plate specification"
    )
    parser.add_argument(
        "-w", "--well-spec", default="B01", help="the well specification"
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        default=None,
        help="the directory containing all output files",
    )
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
                toil,
                options.data_path,
                options.plate_spec,
                [options.well_spec],
                options.source_scheme,
            )

            # Construct and start the NOVOPlasty job
            novoplasty_job = NovoplastyJob(
                read_one_file_ids[0], read_two_file_ids[0], options.output_directory
            )
            novoplasty_rv = toil.start(novoplasty_job)

        else:

            # Restart the NOVOPlasty job
            novoplasty_rv = toil.restart(novoplasty_job)

        # Export the NOVOPlasty log and corrections files, and contigs
        # FASTA file from the file store
        utilities.exportFiles(
            toil, options.output_directory, novoplasty_rv["novoplasty_rv"]
        )
