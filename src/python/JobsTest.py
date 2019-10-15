import filecmp
import os
import shutil
import unittest

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from SpadesJob import SpadesJob
from WellAssemblyJob import WellAssemblyJob
from PlateAssemblyJob import PlateAssemblyJob

import utilities


class ToilTestCase(unittest.TestCase):

    def setUp(self):
        self.data_directory = os.path.join("..", "..", "dat", "miscellaneous")

        self.plate_spec = "A11967A_sW0154"
        self.well_spec = "B01"
        self.coverage_cutoff = "100"

        self.output_directory = "{0}_{1}".format(self.plate_spec, self.well_spec)
        self.actual_directory_a = self.output_directory
        self.test_directory_a = self.output_directory + "_TEST"
        if not os.path.exists(self.test_directory_a):
            os.mkdir(self.test_directory_a)

        self.well_specs = ["G05", "G06"]

        self.actual_directory_b = "A11967B_sW0154"
        self.test_directory_b = self.actual_directory_b + "_TEST"
        if not os.path.exists(self.test_directory_b):
            os.mkdir(self.test_directory_b)

        self.spades_fasta = "contigs.fasta"
        self.apc_fasta = "apc.1.fa"

    def tearDown(self):
        if os.path.exists(self.test_directory_a):
            shutil.rmtree(self.test_directory_a)
        if os.path.exists(self.test_directory_b):
            shutil.rmtree(self.test_directory_b)

    def _import_read_files(self, toil, well_specs):
        read_one_file_ids, read_two_file_ids = utilities.importReadFiles(
            toil, self.data_directory, self.plate_spec, well_specs)
        return read_one_file_ids, read_two_file_ids

    def _import_contigs_file(self, toil):
        contigs_file_id = utilities.importContigsFile(
            toil, self.output_directory)
        return contigs_file_id

    def _assert_true_cmp_spades_fasta(self, test_directory, actual_directory):
        self.assertTrue(
            filecmp.cmp(
                os.path.join(test_directory, self.spades_fasta),
                os.path.join(actual_directory, self.spades_fasta)))

    def _assert_true_cmp_apc_fasta(self, test_directory, actual_directory):
        self.assertTrue(
            filecmp.cmp(
                os.path.join(test_directory, self.apc_fasta),
                os.path.join(actual_directory, self.apc_fasta)))


class JobsTestCase(ToilTestCase):

    def test_spades_job(self):

        options = Job.Runner.getDefaultOptions("spadesJobStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec])

            spades_job = SpadesJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.coverage_cutoff,
                self.output_directory,
            )
            spades_rv = toil.start(spades_job)

            utilities.exportFiles(
                toil, self.test_directory_a, spades_rv['spades_rv'])

        self._assert_true_cmp_spades_fasta(
            self.test_directory_a, self.actual_directory_a)

    def test_apc_job(self):

        options = Job.Runner.getDefaultOptions("apcJobStore")

        with Toil(options) as toil:

            contigs_file_id = self._import_contigs_file(toil)

            apc_job = ApcJob(contigs_file_id)
            apc_rv = toil.start(apc_job)

            utilities.exportFiles(
                toil, self.test_directory_a, apc_rv['apc_rv'])

        self._assert_true_cmp_apc_fasta(
            self.test_directory_a, self.actual_directory_a)


class WellAssemblyJobTestCase(ToilTestCase):

    def test_well_assembly_job(self):

        options = Job.Runner.getDefaultOptions("wellAssemblyJobStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, [self.well_spec])

            well_assembly_job = WellAssemblyJob(
                read_one_file_ids[0],
                read_two_file_ids[0],
                self.coverage_cutoff,
                self.output_directory,
            )
            well_assembly_rv = toil.start(well_assembly_job)

            utilities.exportFiles(
                toil, self.test_directory_a, well_assembly_rv['spades_rv'])
            utilities.exportFiles(
                toil, self.test_directory_a, well_assembly_rv['apc_rv'])

        self._assert_true_cmp_spades_fasta(
            self.test_directory_a, self.actual_directory_a)
        self._assert_true_cmp_apc_fasta(
            self.test_directory_a, self.actual_directory_a)


class PlateAssemblyJobTestCase(ToilTestCase):

    def test_plate_assembly_job(self):

        options = Job.Runner.getDefaultOptions("plateAssemblyJobStore")

        with Toil(options) as toil:

            read_one_file_ids, read_two_file_ids = self._import_read_files(
                toil, self.well_specs)

            plate_assembly_job = PlateAssemblyJob(
                self.well_specs,
                read_one_file_ids,
                read_two_file_ids,
                self.plate_spec,
                self.coverage_cutoff,
                cores=2,
                disk="3G",
                memory="4G",
                )
            well_assembly_rvs = toil.start(plate_assembly_job)

            utilities.exportWellAssemblyFiles(
                toil, self.test_directory_b, self.well_specs, well_assembly_rvs)

        for well_spec in self.well_specs:
            self._assert_true_cmp_spades_fasta(
                os.path.join(self.test_directory_b, well_spec),
                os.path.join(self.actual_directory_b, well_spec))
            self._assert_true_cmp_apc_fasta(
                os.path.join(self.test_directory_b, well_spec),
                os.path.join(self.actual_directory_b, well_spec))


if __name__ == '__main__':

    jobsTestSuite = unittest.TestSuite()
    jobsTestSuite.addTest(JobsTestCase('test_spades_job'))
    jobsTestSuite.addTest(JobsTestCase('test_apc_job'))

    wellJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        WellAssemblyJobTestCase)

    plateJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        PlateJobTestCase)

    testSuites = unittest.TestSuite([
        jobsTestSuite,
        wellJobTestSuite,
        plateJobTestSuite,
    ])

    unittest.TextTestRunner(verbosity=2).run(testSuites)
