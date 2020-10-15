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

import utilities


class ToilTestCase(unittest.TestCase):
    def setUp(self):
        cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-4]
        cmps.extend(["dat", "miscellaneous"])
        self.data_directory = os.sep + os.path.join(*cmps)
        self.plate_spec_a = "A11967A_sW0154"
        self.well_spec_a = "A01"
        self.well_spec_b = "B01"
        self.assembler = "spades"
        self.coverage_cutoff = "100"

        self.output_directory_aa = "{0}_{1}".format(self.plate_spec_a, self.well_spec_a)
        self.actual_directory_aa = self.output_directory_aa
        self.test_directory_aa = self.output_directory_aa + "_TEST"
        if not os.path.exists(self.test_directory_aa):
            os.mkdir(self.test_directory_aa)

        self.output_directory_ab = "{0}_{1}".format(self.plate_spec_a, self.well_spec_b)
        self.actual_directory_ab = self.output_directory_ab
        self.test_directory_ab = self.output_directory_ab + "_TEST"
        if not os.path.exists(self.test_directory_ab):
            os.mkdir(self.test_directory_ab)

        self.plate_spec_b = "A11967B_sW0154"
        self.well_specs = ["G05", "G06"]

        self.actual_directory_bg = self.plate_spec_b
        self.test_directory_bg = self.actual_directory_bg + "_TEST"
        if not os.path.exists(self.test_directory_bg):
            os.mkdir(self.test_directory_bg)

        self.test_masurca_fasta = "final.genome.scf.fasta"
        self.test_novoplasty_fasta = "Circularized_assembly_1_Toil.fasta"
        self.test_shovill_fasta = "contigs.fa"
        self.test_skesa_fasta = "contigs.fa"
        self.test_spades_fasta = "contigs.fasta"
        self.test_unicycler_fasta = "assembly.fasta"

        self.actual_masurca_fasta = "assembly_masurca.fasta"
        self.actual_novoplasty_fasta = "assembly_novoplasty.fasta"
        self.actual_shovill_fasta = "assembly_shovill.fasta"
        self.actual_skesa_fasta = "assembly_skesa.fasta"
        self.actual_spades_fasta = "assembly_spades.fasta"
        self.actual_unicycler_fasta = "assembly_unicycler.fasta"

        self.apc_fasta = "apc.1.fa"

    def tearDown(self):
        # if os.path.exists(self.test_directory_aa):
        #     shutil.rmtree(self.test_directory_aa)
        # if os.path.exists(self.test_directory_ab):
        #     shutil.rmtree(self.test_directory_ab)
        # if os.path.exists(self.test_directory_bg):
        #     shutil.rmtree(self.test_directory_bg)
        pass

    def _import_read_files(self, toil, well_specs):
        read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
            toil, self.data_directory, self.plate_spec_a, well_specs
        )
        return read_one_file_ids, read_two_file_ids

    def _import_contigs_file(self, toil):
        contigs_file_id = utilities.importContigsFile(
            toil, os.path.abspath(self.output_directory_ab)
        )
        return contigs_file_id

    def _assert_true_cmp_fasta(
        self, test_directory, test_fasta, actual_directory, actual_fasta
    ):
        with open(os.path.join(test_directory, test_fasta), "r") as f:
            test_lines = f.readlines()
        with open(os.path.join(actual_directory, actual_fasta), "r") as f:
            actual_lines = f.readlines()
        self.assertTrue(test_lines[1:] == actual_lines[1:])


class JobsTestCase(ToilTestCase):
    def test_masurca_job(self):

        options = Job.Runner.getDefaultOptions("masurcaFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            masurca_job = MasurcaJob(read_one_file_ids[0], read_two_file_ids[0],)
            masurca_rv = toil.start(masurca_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, masurca_rv["masurca_rv"]
            )

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_masurca_fasta,
            self.actual_directory_ab,
            self.actual_masurca_fasta,
        )

    def test_novoplasty_job(self):

        options = Job.Runner.getDefaultOptions("novoplastyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_a]
            )

            novoplasty_job = NovoplastyJob(read_one_file_ids[0], read_two_file_ids[0],)
            novoplasty_rv = toil.start(novoplasty_job)

            utilities.exportFiles(
                toil, self.test_directory_aa, novoplasty_rv["novoplasty_rv"]
            )

        self._assert_true_cmp_fasta(
            self.test_directory_aa,
            self.test_novoplasty_fasta,
            self.actual_directory_aa,
            self.actual_novoplasty_fasta,
        )

    def test_shovill_job(self):

        options = Job.Runner.getDefaultOptions("shovillFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            shovill_job = ShovillJob(
                read_one_file_ids[0], read_two_file_ids[0], self.output_directory_ab,
            )
            shovill_rv = toil.start(shovill_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, shovill_rv["shovill_rv"]
            )

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_shovill_fasta,
            self.actual_directory_ab,
            self.actual_shovill_fasta,
        )

    def test_skesa_job(self):

        options = Job.Runner.getDefaultOptions("skesaFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            skesa_job = SkesaJob(read_one_file_ids[0], read_two_file_ids[0],)
            skesa_rv = toil.start(skesa_job)

            utilities.exportFiles(toil, self.test_directory_ab, skesa_rv["skesa_rv"])

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_skesa_fasta,
            self.actual_directory_ab,
            self.actual_skesa_fasta,
        )

    def test_spades_job(self):

        options = Job.Runner.getDefaultOptions("spadesFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.output_directory_ab,
                {'coverage-cutoff': 100},
            )
            spades_rv = toil.start(spades_job)

            utilities.exportFiles(toil, self.test_directory_ab, spades_rv["spades_rv"])

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_spades_fasta,
            self.actual_directory_ab,
            self.actual_spades_fasta,
        )

    def test_unicycler_job(self):

        options = Job.Runner.getDefaultOptions("unicyclerFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            unicycler_job = UnicyclerJob(
                read_one_file_ids[0], read_two_file_ids[0], self.output_directory_ab,
            )
            unicycler_rv = toil.start(unicycler_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, unicycler_rv["unicycler_rv"]
            )

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_unicycler_fasta,
            self.actual_directory_ab,
            self.actual_unicycler_fasta,
        )

    def test_apc_job(self):

        options = Job.Runner.getDefaultOptions("apcFileStore")

        with Toil(options) as toil:

            contigs_file_id = self._import_contigs_file(toil)

            apc_job = ApcJob(contigs_file_id)
            apc_rv = toil.start(apc_job)

            utilities.exportFiles(toil, self.test_directory_ab, apc_rv["apc_rv"])

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.apc_fasta,
            self.actual_directory_ab,
            self.apc_fasta,
        )


class WellAssemblyJobTestCase(ToilTestCase):
    def test_well_assembly_job(self):

        options = Job.Runner.getDefaultOptions("wellAssemblyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b]
            )

            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.assembler,
                self.coverage_cutoff,
                self.output_directory_ab,
            )
            well_assembly_rv = toil.start(well_assembly_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, well_assembly_rv["spades_rv"]
            )
            utilities.exportFiles(
                toil, self.test_directory_ab, well_assembly_rv["apc_rv"]
            )

        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.test_spades_fasta,
            self.actual_directory_ab,
            self.actual_spades_fasta,
        )
        self._assert_true_cmp_fasta(
            self.test_directory_ab,
            self.apc_fasta,
            self.actual_directory_ab,
            self.apc_fasta,
        )


class PlateAssemblyJobTestCase(ToilTestCase):
    def test_plate_assembly_job(self):

        options = Job.Runner.getDefaultOptions("plateAssemblyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.well_specs
            )

            plate_assembly_job = PlateAssemblyJob(
                self.well_specs,
                read_one_file_ids,
                read_two_file_ids,
                self.plate_spec_b,
                self.assembler,
                self.coverage_cutoff,
                cores=2,
                disk="3G",
                memory="4G",
            )
            well_assembly_rvs = toil.start(plate_assembly_job)

            utilities.exportWellAssemblyFiles(
                toil,
                self.assembler,
                self.test_directory_bg,
                self.well_specs,
                well_assembly_rvs,
            )

        for well_spec in self.well_specs:
            self._assert_true_cmp_fasta(
                os.path.join(self.test_directory_bg, well_spec),
                self.test_spades_fasta,
                os.path.join(self.actual_directory_bg, well_spec),
                self.test_spades_fasta,
            )
            self._assert_true_cmp_fasta(
                os.path.join(self.test_directory_bg, well_spec),
                self.apc_fasta,
                os.path.join(self.actual_directory_bg, well_spec),
                self.apc_fasta,
            )


if __name__ == "__main__":

    jobsTestSuite = unittest.TestSuite()
    # jobsTestSuite.addTest(JobsTestCase("test_masurca_job"))
    # jobsTestSuite.addTest(JobsTestCase("test_novoplasty_job"))
    # jobsTestSuite.addTest(JobsTestCase("test_shovill_job"))
    # jobsTestSuite.addTest(JobsTestCase("test_skesa_job"))
    # jobsTestSuite.addTest(JobsTestCase("test_spades_job"))
    # jobsTestSuite.addTest(JobsTestCase("test_unicycler_job"))
    jobsTestSuite.addTest(JobsTestCase("test_apc_job"))

    # wellAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
    # WellAssemblyJobTestCase
    # )

    # plateAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
    # PlateAssemblyJobTestCase
    # )

    # testSuites = unittest.TestSuite(
    # [jobsTestSuite, wellAssemblyJobTestSuite, plateAssemblyJobTestSuite,]
    # )

    testSuites = unittest.TestSuite([jobsTestSuite,])

    unittest.TextTestRunner(verbosity=2).run(testSuites)
