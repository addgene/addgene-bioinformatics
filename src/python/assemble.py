import os

from toil.job import Job
from toil.common import Toil

from AssemblyJob import AssemblyJob


PLATE_SPEC = "A11967A_sW0154"
WELL_SPEC = "A01"
COVERAGE_CUTOFF = "100"


if __name__ == "__main__":
    parser = Job.Runner.getDefaultArgumentParser()
    options = parser.parse_args()
    options.logLevel = "INFO"

    with Toil(options) as toil:
        if not toil.options.restart:

            # Import the local read files into the file store
            read_one_file_id = toil.importFile("file://" + os.path.abspath(
                "../../dat/miscellaneous/{0}_FASTQ/{0}_{1}_R1_001.fastq.gz".format(
                    PLATE_SPEC, WELL_SPEC)))
            read_two_file_id = toil.importFile("file://" + os.path.abspath(
                "../../dat/miscellaneous/{0}_FASTQ/{0}_{1}_R2_001.fastq.gz".format(
                    PLATE_SPEC, WELL_SPEC)))

            # Construct and start the assembly job
            assembly_job = AssemblyJob(
                read_one_file_id,
                read_two_file_id,
                COVERAGE_CUTOFF,
                PLATE_SPEC,
                )
            assembly_rv = toil.start(assembly_job)

        else:

            # Restart the job
            assembly_rv = toil.restart(assembly_job)

        # Export the SPAdes warnings and log files, and contigs FASTA
        # file from the file store
        spades_rv = assembly_rv['spades_rv']
        toil.exportFile(spades_rv['warnings_file_id'],
                        "file://" + os.path.abspath("./warnings.log"))
        toil.exportFile(spades_rv['spades_file_id'],
                        "file://" + os.path.abspath("./spades.log"))
        toil.exportFile(spades_rv['contigs_file_id'],
                        "file://" + os.path.abspath("./contigs.fasta"))

        # Export the apc output file, and sequence FASTA file from the
        # file store
        apc_rv = assembly_rv['apc_rv']
        toil.exportFile(apc_rv['output_file_id'],
                        "file://" + os.path.abspath("./apc.out"))
        toil.exportFile(apc_rv['sequence_file_id'],
                        "file://" + os.path.abspath("./apc.1.fa"))
