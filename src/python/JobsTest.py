import filecmp
import os
import shutil
import unittest

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from ShovillJob import ShovillJob
from SpadesJob import SpadesJob
from WellAssemblyJob import WellAssemblyJob
from PlateAssemblyJob import PlateAssemblyJob
from NovoplastyJob import NovoplastyJob

import utilities


class ToilTestCase(unittest.TestCase):

    def setUp(self):
        cmps = str(os.path.abspath(__file__)).split(os.sep)[0:-3]
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

        self.novoplasty_fasta = "contigs.fa"
        self.shovill_fasta = "contigs.fa"
        self.spades_fasta = "contigs.fasta"
        self.apc_fasta = "apc.1.fa"

    def tearDown(self):
        if os.path.exists(self.test_directory_aa):
            shutil.rmtree(self.test_directory_aa)
        if os.path.exists(self.test_directory_ab):
            shutil.rmtree(self.test_directory_ab)
        if os.path.exists(self.test_directory_bg):
            shutil.rmtree(self.test_directory_bg)

    def _import_read_files(self, toil, well_specs):
        read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
            toil, self.data_directory, self.plate_spec_a, well_specs)
        return read_one_file_ids, read_two_file_ids

    def _import_contigs_file(self, toil):
        contigs_file_id = utilities.importContigsFile(
            toil, os.path.abspath(self.output_directory_ab))
        return contigs_file_id

    def _assert_true_cmp_fasta(self, test_directory, actual_directory, fasta):
        self.assertTrue(
            filecmp.cmp(
                os.path.join(test_directory, fasta),
                os.path.join(actual_directory, fasta)))


class JobsTestCase(ToilTestCase):

    def test_novoplasty_job(self):
        options = Job.Runner.getDefaultOptions("novoplastyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_a])

            novoplasty_job = NovoplastyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
            )
            novoplasty_rv = toil.start(novoplasty_job)

            utilities.exportFiles(
                toil, self.test_directory_aa, novoplasty_rv['novoplasty_rv'])

        self._assert_true_cmp_fasta(
            self.test_directory_aa, self.actual_directory_aa, self.novoplasty_fasta)

    def test_shovill_job(self):

        options = Job.Runner.getDefaultOptions("shovillFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b])

            shovill_job = ShovillJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.output_directory_ab,
            )
            shovill_rv = toil.start(shovill_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, shovill_rv['shovill_rv'])

        self._assert_true_cmp_fasta(
            self.test_directory_ab, self.actual_directory_ab, self.shovill_fasta)

    def test_spades_job(self):

        options = Job.Runner.getDefaultOptions("spadesFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b])

            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.coverage_cutoff,
                self.output_directory_ab,
            )
            spades_rv = toil.start(spades_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, spades_rv['spades_rv'])

        self._assert_true_cmp_fasta(
            self.test_directory_ab, self.actual_directory_ab, self.spades_fasta)

    def test_apc_job(self):

        options = Job.Runner.getDefaultOptions("apcFileStore")

        with Toil(options) as toil:

            contigs_file_id = self._import_contigs_file(toil)

            apc_job = ApcJob(contigs_file_id)
            apc_rv = toil.start(apc_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, apc_rv['apc_rv'])

        self._assert_true_cmp_fasta(
            self.test_directory_ab, self.actual_directory_ab, self.apc_fasta)


class WellAssemblyJobTestCase(ToilTestCase):

    def test_well_assembly_job(self):

        options = Job.Runner.getDefaultOptions("wellAssemblyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec_b])

            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.assembler,
                self.coverage_cutoff,
                self.output_directory_ab,
            )
            well_assembly_rv = toil.start(well_assembly_job)

            utilities.exportFiles(
                toil, self.test_directory_ab, well_assembly_rv['spades_rv'])
            utilities.exportFiles(
                toil, self.test_directory_ab, well_assembly_rv['apc_rv'])

        self._assert_true_cmp_fasta(
            self.test_directory_ab, self.actual_directory_ab, self.spades_fasta)
        self._assert_true_cmp_fasta(
            self.test_directory_ab, self.actual_directory_ab, self.apc_fasta)


class PlateAssemblyJobTestCase(ToilTestCase):

    def test_plate_assembly_job(self):

        options = Job.Runner.getDefaultOptions("plateAssemblyFileStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.well_specs)

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
                toil, self.assembler, self.test_directory_bg, self.well_specs, well_assembly_rvs)

        for well_spec in self.well_specs:
            self._assert_true_cmp_fasta(
                os.path.join(self.test_directory_bg, well_spec),
                os.path.join(self.actual_directory_bg, well_spec),
                self.spades_fasta)
            self._assert_true_cmp_fasta(
                os.path.join(self.test_directory_bg, well_spec),
                os.path.join(self.actual_directory_bg, well_spec),
                self.apc_fasta)


if __name__ == '__main__':

    jobsTestSuite = unittest.TestSuite()
    jobsTestSuite.addTest(JobsTestCase("test_novoplasty_job"))
    jobsTestSuite.addTest(JobsTestCase("test_shovill_job"))
    jobsTestSuite.addTest(JobsTestCase("test_spades_job"))
    jobsTestSuite.addTest(JobsTestCase("test_apc_job"))

    wellAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        WellAssemblyJobTestCase)

    plateAssemblyJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        PlateAssemblyJobTestCase)

    testSuites = unittest.TestSuite([
        jobsTestSuite,
        wellAssemblyJobTestSuite,
        plateAssemblyJobTestSuite,
    ])

    unittest.TextTestRunner(verbosity=2).run(testSuites)
