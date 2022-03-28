import os
import shutil
import unittest

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from MasurcaJob import MasurcaJob
from NovoplastyJob import NovoplastyJob
from PlateAssemblyJob import PlateAssemblyJob
from ShovillJob import ShovillJob
from SpadesJob import SpadesJob
from SkesaJob import SkesaJob
from UnicyclerJob import UnicyclerJob
from WellAssemblyJob import WellAssemblyJob
from BBDukJob import BBDukJob
from BBNormJob import BBNormJob
from BBMergeJob import BBMergeJob

import utilities


class ToilTestCase(unittest.TestCase):
    def setUp(self):

        cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
        cmps.extend(["dat", "miscellaneous"])
        self.data_path = os.sep + os.path.join(*cmps)

        self.plate_spec_a = "A11967A_sW0154"
        self.plate_spec_b = "A11967B_sW0154"
        self.plate_spec_c = "A14735A_sW0248"

        self.well_spec_a = "A01"
        self.well_spec_b = "B01"
        self.well_specs = ["G04", "G06"]

        cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-1]
        self.config_path = os.sep + os.path.join(*cmps)
        self.config_file = "Assembler.ini"
        self.adapters_path = os.sep + os.path.join(*cmps)
        self.adapters_file = "adapters.fa"

        self.test_directory = os.path.join(os.sep + os.path.join(*cmps), "test")
        if not os.path.exists(self.test_directory):
            os.mkdir(self.test_directory)
        self.output_directory = os.path.join(os.sep + os.path.join(*cmps), "output")

    def tearDown(self):
        if os.path.exists(self.test_directory):
            shutil.rmtree(self.test_directory)

    def _import_read_files(self, toil, plate_spec, well_specs):
        read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
            toil, self.data_path, plate_spec, well_specs
        )
        return read_one_file_ids, read_two_file_ids

    def _import_config_file(self, toil):
        config_file_id = utilities.importConfigFile(
            toil, os.path.join(self.config_path, self.config_file)
        )
        return config_file_id

    def _import_adapters_file(self, toil):
        adapters_file_id = utilities.importAdaptersFile(
            toil, os.path.join(self.adapters_path, self.adapters_file)
        )
        return adapters_file_id

    def _assert_true_cmp_fasta(
        self, test_directory, test_fasta, actual_directory, actual_fasta
    ):
        with open(os.path.join(test_directory, test_fasta), "r") as f:
            test_lines = f.readlines()
        with open(os.path.join(actual_directory, actual_fasta), "r") as f:
            actual_lines = f.readlines()
        self.assertTrue(test_lines[1:] == actual_lines[1:])

    def _assert_true_cmp_sorted_fasta(
        self, test_directory, test_fasta, actual_directory, actual_fasta
    ):
        with open(os.path.join(test_directory, test_fasta), "r") as f:
            test_lines = f.readlines()
        with open(os.path.join(actual_directory, actual_fasta), "r") as f:
            actual_lines = f.readlines()
        test_lines.sort()
        actual_lines.sort()
        self.assertTrue(test_lines[1:] == actual_lines[1:])


class JobsTestCase(ToilTestCase):
    def test_masurca_job(self):

        options = Job.Runner.getDefaultOptions("masurcaFileStore")

        test_directory = os.path.join(
            self.test_directory, "MasurcaJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "MasurcaJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)

            masurca_job = MasurcaJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
            )
            masurca_rv = toil.start(masurca_job)

            utilities.exportFiles(toil, test_directory, masurca_rv["masurca_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "final.genome.scf.fasta",
            output_directory,
            "final.genome.scf.fasta",
        )

    def test_novoplasty_job(self):

        options = Job.Runner.getDefaultOptions("novoplastyFileStore")

        test_directory = os.path.join(
            self.test_directory, "NovoplastyJob", self.plate_spec_a, self.well_spec_a
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "NovoplastyJob", self.plate_spec_a, self.well_spec_a
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_a]
            )

            config_file_id = self._import_config_file(toil)

            novoplasty_job = NovoplastyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
            )
            novoplasty_rv = toil.start(novoplasty_job)

            utilities.exportFiles(toil, test_directory, novoplasty_rv["novoplasty_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "Circularized_assembly_1_Toil.fasta",
            output_directory,
            "Circularized_assembly_1_Toil.fasta",
        )

    def test_shovill_job(self):

        options = Job.Runner.getDefaultOptions("shovillFileStore")

        test_directory = os.path.join(
            self.test_directory, "ShovillJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "ShovillJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)

            shovill_job = ShovillJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
                os.path.join("test", "ShovillJob", self.plate_spec_a, self.well_spec_b),
            )
            shovill_rv = toil.start(shovill_job)

            utilities.exportFiles(toil, test_directory, shovill_rv["shovill_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "contigs.fa",
            output_directory,
            "contigs.fa",
        )

    def test_skesa_job(self):

        options = Job.Runner.getDefaultOptions("skesaFileStore")

        test_directory = os.path.join(
            self.test_directory, "SkesaJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "SkesaJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)

            skesa_job = SkesaJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
            )
            skesa_rv = toil.start(skesa_job)

            utilities.exportFiles(toil, test_directory, skesa_rv["skesa_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "contigs.fa",
            output_directory,
            "contigs.fa",
        )

    def test_spades_job(self):

        options = Job.Runner.getDefaultOptions("spadesFileStore")

        test_directory = os.path.join(
            self.test_directory, "SpadesJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "SpadesJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)

            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
                os.path.join("test", "SpadesJob", self.plate_spec_a, self.well_spec_b),
            )
            spades_rv = toil.start(spades_job)

            utilities.exportFiles(toil, test_directory, spades_rv["spades_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "contigs.fasta",
            output_directory,
            "contigs.fasta",
        )

    def test_unicycler_job(self):

        options = Job.Runner.getDefaultOptions("unicyclerFileStore")

        test_directory = os.path.join(
            self.test_directory, "UnicyclerJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "UnicyclerJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)

            unicycler_job = UnicyclerJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
                os.path.join(
                    "test", "UnicyclerJob", self.plate_spec_a, self.well_spec_b
                ),
            )
            unicycler_rv = toil.start(unicycler_job)

            utilities.exportFiles(toil, test_directory, unicycler_rv["unicycler_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "assembly.fasta",
            output_directory,
            "assembly.fasta",
        )

    def test_apc_job(self):

        options = Job.Runner.getDefaultOptions("apcFileStore")

        test_directory = os.path.join(
            self.test_directory, "ApcJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "ApcJob", self.plate_spec_a, self.well_spec_b
        )

        with Toil(options) as toil:

            contigs_file_id = utilities.importFile(
                toil,
                os.path.join(
                    self.output_directory,
                    "SpadesJob",
                    self.plate_spec_a,
                    self.well_spec_b,
                    "contigs.fasta",
                ),
            )

            apc_job = ApcJob(contigs_file_id)
            apc_rv = toil.start(apc_job)

            utilities.exportFiles(toil, test_directory, apc_rv["apc_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "apc.1.fa",
            output_directory,
            "apc.1.fa",
        )

    def test_bbduk_job(self):

        options = Job.Runner.getDefaultOptions("bbdukFileStore")

        test_directory = os.path.join(
            self.test_directory, "BBDukJob", self.plate_spec_c, self.well_spec_a
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "BBDukJob", self.plate_spec_c, self.well_spec_a
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_c, [self.well_spec_a]
            )

            config_file_id = self._import_config_file(toil)
            adapters_file_id = self._import_adapters_file(toil)

            bbduk_job = BBDukJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
                adapters_file_id,
                self.adapters_file,
            )
            bbduk_rv = toil.start(bbduk_job)

            utilities.exportFiles(toil, test_directory, bbduk_rv["bbduk_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "trimmed_R1.fastq",
            output_directory,
            "trimmed_R1.fastq",
        )

        self._assert_true_cmp_fasta(
            test_directory,
            "trimmed_R2.fastq",
            output_directory,
            "trimmed_R2.fastq",
        )

    def test_bbnorm_job(self):

        options = Job.Runner.getDefaultOptions("bbnormFileStore")

        test_directory = os.path.join(
            self.test_directory, "BBNormJob", self.plate_spec_c, self.well_spec_a
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.test_directory, "BBNormJob", self.plate_spec_c, self.well_spec_a
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_c, [self.well_spec_a]
            )

            config_file_id = self._import_config_file(toil)

            bbnorm_job = BBNormJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
            )
            bbnorm_rv = toil.start(bbnorm_job)

            utilities.exportFiles(toil, test_directory, bbnorm_rv["bbnorm_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "normed_R1.fastq",
            output_directory,
            "normed_R1.fastq",
        )

        self._assert_true_cmp_fasta(
            test_directory,
            "normed_R2.fastq",
            output_directory,
            "normed_R2.fastq",
        )

    def test_bbmerge_job(self):

        options = Job.Runner.getDefaultOptions("bbmergeFileStore")

        test_directory = os.path.join(
            self.test_directory, "BBMerge", self.plate_spec_c, self.well_spec_a
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory, "BBMergeJob", self.plate_spec_c, self.well_spec_a
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_c, [self.well_spec_a]
            )

            config_file_id = self._import_config_file(toil)

            bbmerge_job = BBMergeJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                config_file_id,
                self.config_file,
            )
            bbmerge_rv = toil.start(bbmerge_job)

            utilities.exportFiles(toil, test_directory, bbmerge_rv["bbmerge_rv"])

        self._assert_true_cmp_sorted_fasta(
            test_directory,
            "unmerged_R1.fastq",
            output_directory,
            "unmerged_R1.fastq",
        )

        self._assert_true_cmp_sorted_fasta(
            test_directory,
            "unmerged_R2.fastq",
            output_directory,
            "unmerged_R2.fastq",
        )

        self._assert_true_cmp_sorted_fasta(
            test_directory,
            "merged.fastq",
            output_directory,
            "merged.fastq",
        )


class WellAssemblyJobTestCase(ToilTestCase):
    def test_well_assembly_job(self):

        options = Job.Runner.getDefaultOptions("wellAssemblyFileStore")

        test_directory = os.path.join(
            self.test_directory, "WellAssemblyJob", self.plate_spec_a, self.well_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)
        output_directory = os.path.join(
            self.output_directory,
            "WellAssemblyJob",
            self.plate_spec_a,
            self.well_spec_b,
        )

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_a, [self.well_spec_b]
            )

            config_file_id = self._import_config_file(toil)
            adapters_file_id = self._import_adapters_file(toil)

            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                "spades",
                config_file_id,
                self.config_file,
                adapters_file_id,
                self.adapters_file,
                os.path.join(
                    "test", "WellAssemblyJob", self.plate_spec_a, self.well_spec_b
                ),
                preprocessing=False,
                export_intermediary_files=True,
            )
            well_assembly_rv = toil.start(well_assembly_job)

            utilities.exportFiles(toil, test_directory, well_assembly_rv["spades_rv"])
            utilities.exportFiles(toil, test_directory, well_assembly_rv["apc_rv"])

        self._assert_true_cmp_fasta(
            test_directory,
            "contigs.fasta",
            output_directory,
            "contigs.fasta",
        )
        self._assert_true_cmp_fasta(
            test_directory,
            "apc.1.fa",
            output_directory,
            "apc.1.fa",
        )


class PlateAssemblyJobTestCase(ToilTestCase):
    def test_plate_assembly_job(self):

        options = Job.Runner.getDefaultOptions("plateAssemblyFileStore")

        test_directory = os.path.join(
            self.test_directory, "PlateAssemblyJob", self.plate_spec_b
        )
        os.makedirs(test_directory, exist_ok=True)

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.plate_spec_b, self.well_specs
            )

            config_file_id = self._import_config_file(toil)
            adapters_file_id = self._import_adapters_file(toil)

            plate_assembly_job = PlateAssemblyJob(
                self.well_specs,
                read_one_file_ids,
                read_two_file_ids,
                self.plate_spec_a,
                "spades",
                config_file_id,
                self.config_file,
                adapters_file_id,
                self.adapters_file,
                preprocessing=False,
                cores=2,
                disk="3G",
                memory="4G",
            )
            well_assembly_rvs = toil.start(plate_assembly_job)

            utilities.exportWellAssemblyFiles(
                toil,
                "spades",
                test_directory,
                self.well_specs,
                well_assembly_rvs,
            )

        for well_spec in self.well_specs:
            test_directory = os.path.join(
                self.test_directory, "PlateAssemblyJob", self.plate_spec_b, well_spec
            )
            os.makedirs(test_directory, exist_ok=True)
            output_directory = os.path.join(
                self.output_directory, "PlateAssemblyJob", self.plate_spec_b, well_spec
            )

            self._assert_true_cmp_fasta(
                test_directory,
                "contigs.fasta",
                output_directory,
                "contigs.fasta",
            )
            self._assert_true_cmp_fasta(
                test_directory,
                "apc.1.fa",
                output_directory,
                "apc.1.fa",
            )


if __name__ == "__main__":

    jobsTestSuite = unittest.TestSuite()

    # TODO: Resolve
    jobsTestSuite.addTest(JobsTestCase("test_masurca_job"))
    jobsTestSuite.addTest(JobsTestCase("test_novoplasty_job"))
    jobsTestSuite.addTest(JobsTestCase("test_shovill_job"))
    jobsTestSuite.addTest(JobsTestCase("test_skesa_job"))
    jobsTestSuite.addTest(JobsTestCase("test_spades_job"))
    jobsTestSuite.addTest(JobsTestCase("test_unicycler_job"))
    jobsTestSuite.addTest(JobsTestCase("test_apc_job"))
    jobsTestSuite.addTest(JobsTestCase("test_bbduk_job"))
    jobsTestSuite.addTest(JobsTestCase("test_bbnorm_job"))
    jobsTestSuite.addTest(JobsTestCase("test_bbmerge_job"))

    wellAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        WellAssemblyJobTestCase
    )

    plateAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        PlateAssemblyJobTestCase
    )

    testSuites = unittest.TestSuite(
        [
            jobsTestSuite,
            wellAssemblyJobTestSuite,
            plateAssemblyJobTestSuite,
        ]
    )

    unittest.TextTestRunner(verbosity=2).run(testSuites)
