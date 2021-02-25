from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil
from toil.lib.docker import apiDockerCall

import utilities

logger = logging.getLogger(__name__)


class MasurcaJob(Job):
    """
    Accepts paired-end Illumina reads for assembly using MaSuRCA.
    """

    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        config_file_id,
        config_file_name,
        parent_rv={},
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
        config_file_id : toil.fileStore.FileID
            id of the file in the file store containing assembler args
        config_file_name : str
            name of the file in the file store containing assembler args
        parent_rv : dict
            dictionary of return values from the parent job
        """
        super(MasurcaJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name
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
        contigs_file_name = "final.genome.scf.fasta"  # Final assembly scaffolds file

        try:
            # Read the config file from the file store into the local
            # temporary directory, and parse
            config_file_path = utilities.readGlobalFile(
                fileStore, self.config_file_id, self.config_file_name
            )
            common_config, assembler_params = utilities.parseConfigFile(
                config_file_path, "masurca"
            )

            # Read the read files from the file store into the local
            # temporary directory
            read_one_file_path = utilities.readGlobalFile(
                fileStore, self.read_one_file_id, common_config["read_one_file_name"]
            )
            read_two_file_path = utilities.readGlobalFile(
                fileStore, self.read_two_file_id, common_config["read_two_file_name"]
            )

            # Write the MaSuRCA config file into the local temporary
            # directory
            working_dir = fileStore.localTempDir
            logger.info("Handling configuration file {0}".format(assembler_params['config_file_name']))
            with open(os.path.join(working_dir, assembler_params['config_file_name']), "w+") as f:
                config = """
# Input file for 'masurca' command to create 'assemble.sh'.
# --------------------------------------------------------------------------------
# DATA is specified as type (PE, JUMP, OTHER, PACBIO) and 5 fields:
#   1) two_letter_prefix
#   2) mean
#   3) stdev
#   4) fastq(.gz)_fwd_reads
#   5) fastq(.gz)_rev_reads
# The PE reads are always assumed to be innies, i.e. --->.<---, and
# JUMP are assumed to be outties <---.--->. If there are any jump
# libraries that are innies, such as longjump, specify them as JUMP
# and specify NEGATIVE mean. Reverse reads are optional for PE
# libraries and mandatory for JUMP libraries. Any OTHER sequence data
# (454, Sanger, Ion torrent, etc) must be first converted into Celera
# Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
# Illumina paired end reads supplied as:
#   <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
# If single-end, do not specify <reverse_reads>
# MUST HAVE Illumina paired end reads to use MaSuRCA
PE = pe 500 50 {read_one_file_path} {read_two_file_path}
#
# Illumina mate pair reads supplied as:
#   <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#JUMP = sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#
# Pacbio OR Nanopore reads must be in a single fasta or fastq file
# with absolute path, can be gzipped. If you have both types of reads
# supply them both as NANOPORE type.
#PACBIO = /FULL_PATH/pacbio.fa
#NANOPORE = /FULL_PATH/nanopore.fa
#
# Other reads (Sanger, 454, etc) one frg file, concatenate your frg
# files into one if you have many
#OTHER = /FULL_PATH/file.frg
END
# --------------------------------------------------------------------------------
PARAMETERS
# Set this to 1 if your Illumina jumping library reads are shorter
# than 100bp
EXTEND_JUMP_READS = 0
# This is k-mer size for deBruijn graph values between 25 and 127 are
# supported, auto will compute the optimal size based on the read data
# and GC content
GRAPH_KMER_SIZE = auto
# Set this to 1 for all Illumina-only assemblies, or set this to 0 if
# you have more than 15x coverage by long reads (Pacbio or Nanopore)
# or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc).
USE_LINKING_MATES = 1
# Specifies whether to run mega-reads correction on the grid
USE_GRID = 0
# Specifies grid engine to use SGE or SLURM
GRID_ENGINE = SGE
# Specifies queue (for SGE) or partition (for SLURM) to use when
# running on the grid MANDATORY
GRID_QUEUE = all.q
# Batch size in the amount of long read sequence for each batch on the
# grid
GRID_BATCH_SIZE = 300000000
# Use at most this much coverage by the longest Pacbio or Nanopore
# reads, discard the rest of the reads
LHE_COVERAGE = 25
# Set to 1 to only do one pass of mega-reads, for faster but worse
# quality assembly
MEGA_READS_ONE_PASS = 0
# This parameter is useful if you have too many Illumina jumping
# library mates. Typically set it to 60 for bacteria and 300 for the
# other organisms
LIMIT_JUMP_COVERAGE = 300
# These are the additional parameters to Celera Assembler. Do not
# worry about performance, number or processors or batch sizes --
# these are computed automatically. Set cgwErrorRate=0.25 for bacteria
# and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS = cgwErrorRate=0.15
# Minimum count k-mers used in error correction 1 means all k-mers are
# used. One can increase to 2 if Illumina coverage > 100
KMER_COUNT_THRESHOLD = 1
# Whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS = 1
# Auto-detected number of cpus to use
NUM_THREADS = {threads}
# This is mandatory jellyfish hash size -- a safe value is
# estimated_genome_size*estimated_coverage
JF_SIZE = 2000000
# Set this to 1 to use SOAPdenovo contigging/scaffolding
# module. Assembly will be worse but will run faster. Useful for very
# large (> 5Gbp) genomes from Illumina-only data.
SOAP_ASSEMBLY = 0
END
# --------------------------------------------------------------------------------
                 """.format(
                    read_one_file_path=read_one_file_path,
                    read_two_file_path=read_two_file_path,
                    threads=assembler_params['threads'],
                )
                f.write(config)

            # Mount the Toil local temporary directory to the same path in
            # the container, and use the path as the working directory in
            # the container, then call MaSuRCA
            # TODO: Specify the container on construction
            image = "ralatsdio/masurca:v3.4.1"
            logger.info("Calling image {0}".format(image))
            parameters = ["masurca.sh", assembler_params['config_file_name'],],
            logger.info("Using parameters {0}".format(str(parameters)))
            logger.info(config)
            import pdb; pdb.set_trace()
            apiDockerCall(
                self,
                image=image,
                volumes={working_dir: {"bind": working_dir, "mode": "rw"}},
                working_dir=working_dir,
                parameters=parameters,
            )

            # Write the CABOG stdout, super read code stderr, and
            # final assembly scaffolds FASTA files from the local
            # temporary directory into the file store
            ca0_file_id = utilities.writeGlobalFile(fileStore, ca0_file_name)
            ca1_file_id = utilities.writeGlobalFile(fileStore, ca1_file_name)
            ca2_file_id = utilities.writeGlobalFile(fileStore, ca2_file_name)
            sr1_file_id = utilities.writeGlobalFile(fileStore, sr1_file_name)
            contigs_file_id = utilities.writeGlobalFile(
                fileStore, "CA", contigs_file_name
            )

        except Exception as exc:
            # Ensure expectred return values on exceptions
            logger.info("Calling image {0} failed: {1}".format(image, exc))
            ca0_file_id = None
            ca1_file_id = None
            ca2_file_id = None
            sr1_file_id = None
            contigs_file_id = None

        # Return file ids and names for export
        masurca_rv = {
            "masurca_rv": {
                "ca0_file": {"id": ca0_file_id, "name": ca0_file_name,},
                "ca1_file": {"id": ca1_file_id, "name": ca1_file_name,},
                "ca2_file": {"id": ca2_file_id, "name": ca2_file_name,},
                "sr1_file": {"id": sr1_file_id, "name": sr1_file_name,},
                "contigs_file": {"id": contigs_file_id, "name": contigs_file_name,},
            }
        }
        masurca_rv.update(self.parent_rv)
        logger.info("Return value {0}".format(masurca_rv))
        return masurca_rv


if __name__ == "__main__":
    """
    Assemble reads corresponding to a single well.
    """
    # Parse FASTQ data path, plate and well specification,
    # configuration path and file, and output directory, making the
    # output directory if needed
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
        "-l", "--plate-spec", default="A11967A_sW0154", help="the plate specification"
    )
    parser.add_argument(
        "-w", "--well-spec", default="B01", help="the well specification"
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

            # Import local config file into the file store
            config_file_id = utilities.importConfigFile(
                toil, os.path.join(options.config_path, options.config_file)
            )

            # Construct and start the MaSuRCA job
            masurca_job = MasurcaJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                options.config_file,
            )
            masurca_rv = toil.start(masurca_job)

        else:

            # Restart the MaSuRCA job
            masurca_rv = toil.restart(masurca_job)

        # Export all MaSuRCA output files from the file store
        utilities.exportFiles(toil, options.output_directory, masurca_rv["masurca_rv"])
