from argparse import ArgumentParser
import logging
import os

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from MasurcaJob import MasurcaJob
from NovoplastyJob import NovoplastyJob
from SpadesJob import SpadesJob
from ShovillJob import ShovillJob
from SkesaJob import SkesaJob
from UnicyclerJob import UnicyclerJob
from BBDukJob import BBDukJob
from BBNormJob import BBNormJob
from BBMergeJob import BBMergeJob

import utilities

logger = logging.getLogger(__name__)


class WellAssemblyJob(Job):
    """
    Runs a SpadesJob or ShovillJob followed by an ApcJob, or a
    NovoplastyJob
    """

    def __init__(
        self,
        read_one_file_id,
        read_two_file_id,
        assembler,
        config_file_id,
        config_file_name,
        adapters_file_id,
        adapters_file_name,
        output_directory,
        preprocessing=True,
        *args,
        **kwargs
    ):
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
            name of assembler to run, from utilities.ASSEMBLERS_TO_RUN
        config_file_id : toil.fileStore.FileID
            id of the file in the file store containing assembler args
        config_file_name : str
            name of the file in the file store containing assembler args
        output_directory : str
            name of directory for output
        """
        super(WellAssemblyJob, self).__init__(*args, **kwargs)
        self.read_one_file_id = read_one_file_id
        self.read_two_file_id = read_two_file_id
        if assembler not in utilities.ASSEMBLERS_TO_RUN:
            raise Exception("Unexpected assembler")
        self.assembler = assembler
        self.config_file_id = config_file_id
        self.config_file_name = config_file_name
        self.adapters_file_id = adapters_file_id
        self.adapters_file_name = adapters_file_name
        self.output_directory = output_directory
        self.preprocessing = preprocessing

    def run(self, fileStore):
        """
        Returns
        -------
        dict
            the return values of the SpadesJob and ApcJob
        """
        try:
            if self.assembler == "masurca":
                masurca_job = MasurcaJob(
                    self.read_one_file_id,
                    self.read_two_file_id,
                    self.config_file_id,
                    self.config_file_name,
                )
                apc_job = ApcJob(
                    masurca_job.rv("masurca_rv", "contigs_file", "id"),
                    parent_rv=masurca_job.rv(),
                )
                final_job = self.addChild(masurca_job).addChild(apc_job)

            elif self.assembler == "novoplasty":
                novoplasty_job = NovoplastyJob(
                    self.read_one_file_id,
                    self.read_two_file_id,
                    self.config_file_id,
                    self.config_file_name,
                )
                final_job = self.addChild(novoplasty_job)

            elif self.assembler == "shovill":
                shovill_job = ShovillJob(
                    self.read_one_file_id,
                    self.read_two_file_id,
                    self.config_file_id,
                    self.config_file_name,
                    self.output_directory,
                )
                apc_job = ApcJob(
                    shovill_job.rv("shovill_rv", "contigs_file", "id"),
                    parent_rv=shovill_job.rv(),
                )
                final_job = self.addChild(shovill_job).addChild(apc_job)

            elif self.assembler == "skesa":
                skesa_job = SkesaJob(
                    self.read_one_file_id,
                    self.read_two_file_id,
                    self.config_file_id,
                    self.config_file_name,
                )
                apc_job = ApcJob(
                    skesa_job.rv("skesa_rv", "contigs_file", "id"),
                    parent_rv=skesa_job.rv(),
                )
                final_job = self.addChild(skesa_job).addChild(apc_job)

            elif self.assembler == "spades":

                if self.preprocessing:
                    ## BBTools preprocessing
                    bbduk_job = BBDukJob(
                        self.read_one_file_id,
                        self.read_two_file_id,
                        self.config_file_id,
                        self.config_file_name,
                        self.adapters_file_id,
                        self.adapters_file_name,
                    )
                    bbnorm_job = BBNormJob(
                        bbduk_job.rv("bbduk_rv", "out1_file", "id"),
                        bbduk_job.rv("bbduk_rv", "out2_file", "id"),
                        self.config_file_id,
                        self.config_file_name,
                        chained_job=True,
                        parent_rv=bbduk_job.rv(),
                    )
                    bbmerge_job = BBMergeJob(
                        bbnorm_job.rv("bbnorm_rv", "out1_file", "id"),
                        bbnorm_job.rv("bbnorm_rv", "out2_file", "id"),
                        self.config_file_id,
                        self.config_file_name,
                        chained_job=True,
                        parent_rv=bbnorm_job.rv(),
                    )
                    spades_job = SpadesJob(
                        bbmerge_job.rv("bbmerge_rv", "outu1_file", "id"),
                        bbmerge_job.rv("bbmerge_rv", "outu2_file", "id"),
                        self.config_file_id,
                        self.config_file_name,
                        self.output_directory,
                        merged_file_id=bbmerge_job.rv(
                            "bbmerge_rv", "merged_file", "id"
                        ),
                        chained_job=True,
                        parent_rv=bbmerge_job.rv(),
                    )
                    apc_job = ApcJob(
                        spades_job.rv("spades_rv", "contigs_file", "id"),
                        parent_rv=spades_job.rv(),
                    )

                    final_job = (
                        self.addChild(bbduk_job)
                        .addChild(bbnorm_job)
                        .addChild(bbmerge_job)
                        .addChild(spades_job)
                        .addChild(apc_job)
                    )

                else:

                    spades_job = SpadesJob(
                        self.read_one_file_id,
                        self.read_two_file_id,
                        self.config_file_id,
                        self.config_file_name,
                        self.output_directory,
                    )
                    apc_job = ApcJob(
                        spades_job.rv("spades_rv", "contigs_file", "id"),
                        parent_rv=spades_job.rv(),
                    )

                    final_job = self.addChild(spades_job).addChild(apc_job)

            elif self.assembler == "unicycler":

                if self.preprocessing:
                    ## BBTools preprocessing
                    bbduk_job = BBDukJob(
                        self.read_one_file_id,
                        self.read_two_file_id,
                        self.config_file_id,
                        self.config_file_name,
                        self.adapters_file_id,
                        self.adapters_file_name,
                    )
                    bbnorm_job = BBNormJob(
                        bbduk_job.rv("bbduk_rv", "out1_file", "id"),
                        bbduk_job.rv("bbduk_rv", "out2_file", "id"),
                        self.config_file_id,
                        self.config_file_name,
                        chained_job=True,
                        parent_rv=bbduk_job.rv(),
                    )
                    unicycler_job = UnicyclerJob(
                        bbnorm_job.rv("bbnorm_rv", "out1_file", "id"),
                        bbnorm_job.rv("bbnorm_rv", "out2_file", "id"),
                        self.config_file_id,
                        self.config_file_name,
                        self.output_directory,
                    )
                    apc_job = ApcJob(
                        unicycler_job.rv("unicycler_rv", "contigs_file", "id"),
                        parent_rv=unicycler_job.rv(),
                    )
                    final_job = (
                        self.addChild(bbduk_job)
                        .addChild(bbnorm_job)
                        .addChild(unicycler_job)
                        .addChild(apc_job)
                    )
                else:
                    unicycler_job = UnicyclerJob(
                        self.read_one_file_id,
                        self.read_two_file_id,
                        self.config_file_id,
                        self.config_file_name,
                        self.output_directory,
                    )
                    apc_job = ApcJob(
                        unicycler_job.rv("unicycler_rv", "contigs_file", "id"),
                        parent_rv=unicycler_job.rv(),
                    )
                    final_job = self.addChild(unicycler_job).addChild(apc_job)

            # Assign assembler return values
            assembler_rv = final_job.rv()

        except Exception as exc:
            # Ensure expected return values on exceptions
            logger.info(
                "Assembler {0} failed for {1}: {2}".format(
                    self.assembler, self.output_directory, exc
                )
            )
            assembler_rv = {self.assembler: {}}

        # Return file ids, and names for export
        return final_job.rv()


if __name__ == "__main__":
    """
    Assemble reads and circularize contigs corresponding to a single
    well.
    """
    # Parse FASTQ data path, plate and well specification, assembler,
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
    parser.add_argument(
        "-a",
        "--assembler",
        default="spades",
        choices=utilities.ASSEMBLERS_TO_RUN,
        help="name of assembler to run",
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
        "--adapters-path",
        default=os.sep + os.path.join(*cmps),
        help="path to a .fa file with adapters to be passed to bbtools",
    )
    parser.add_argument(
        "--adapters-file",
        default="adapters.fa",
        help="path to a .fa file with adapters to be passed to bbtools",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        default=None,
        help="the directory containing all output files",
    )
    parser.add_argument(
        "--no-preprocessing", dest="preprocessing", action="store_false"
    )
    parser.set_defaults(preprocessing=True)

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

            # Import local adapters file into the file store
            adapters_file_id = utilities.importAdaptersFile(
                toil, os.path.join(options.adapters_path, options.adapters_file)
            )

            # Construct and start the well assembly job
            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                options.assembler,
                config_file_id,
                options.config_file,
                adapters_file_id,
                options.adapters_file,
                options.output_directory,
                preprocessing=options.preprocessing,
            )
            well_assembly_rv = toil.start(well_assembly_job)

        else:

            # Restart the well assembly job
            well_assembly_rv = toil.restart(well_assembly_job)

        # Export all assembler output files, and apc output files, if
        # apc is required by the asembler, from the file store
        utilities.exportWellAssemblyFiles(
            toil,
            options.assembler,
            options.plate_spec,
            [options.well_spec],
            [well_assembly_rv],
        )
