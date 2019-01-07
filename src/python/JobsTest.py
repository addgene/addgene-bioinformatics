import filecmp
import os
import unittest

from toil.job import Job
from toil.common import Toil

from ApcJob import ApcJob
from SpadesJob import SpadesJob

import utilities


class JobsTestCase(unittest.TestCase):

    def test_spades_job(self):

        options = Job.Runner.getDefaultOptions("spadesJobStore")

        with Toil(options) as toil:

            read_one_file_id = utilities.importFile(
                toil, os.path.join(
                    "..", "..", "dat", "miscellaneous", "A11967A_sW0154_FASTQ",
                    "A11967A_sW0154_A01_R1_001.fastq.gz"
                )
            )
            read_two_file_id = utilities.importFile(
                toil, os.path.join(
                    "..", "..", "dat", "miscellaneous", "A11967A_sW0154_FASTQ",
                    "A11967A_sW0154_A01_R2_001.fastq.gz"
                )
            )

            coverage_cutoff = "100"
            output_directory = "A11967A_sW0154_A01"

            spades_job = SpadesJob(
                read_one_file_id,
                read_two_file_id,
                coverage_cutoff,
                output_directory,
            )
            spades_rv = toil.start(spades_job)

            utilities.exportFiles(toil, spades_rv['spades_rv'])

        self.assertTrue(
            filecmp.cmp(
                os.path.join("A11967A_sW0154_A01", "contigs.fasta"),
                "contigs.fasta"
            )
        )

    def test_apc_job(self):

        options = Job.Runner.getDefaultOptions("apcJobStore")

        with Toil(options) as toil:

            contigs_file_id = utilities.importFile(
                toil, os.path.join("A11967A_sW0154_A01", "contigs.fasta"))

            apc_job = ApcJob(contigs_file_id)
            apc_rv = toil.start(apc_job)

            utilities.exportFiles(toil, apc_rv['apc_rv'])

        self.assertTrue(
            filecmp.cmp(
                os.path.join("A11967A_sW0154_A01", "apc.1.fa"),
                "apc.1.fa"
            )
        )


class WellJobTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_well_assembly_job(self):
        pass


class PlateJobTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_plate_assembly_job(self):
        pass


if __name__ == '__main__':

    jobsTestSuite = unittest.TestSuite()
    jobsTestSuite.addTest(JobsTestCase('test_spades_job'))
    jobsTestSuite.addTest(JobsTestCase('test_apc_job'))

    wellJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        WellJobTestCase)

    plateJobTestSuite = unittest.TestLoader().loadTestsFromTestCase(
        PlateJobTestCase)

    testSuites = unittest.TestSuite([
        jobsTestSuite,
        wellJobTestSuite,
        plateJobTestSuite,
    ])

    unittest.TextTestRunner(verbosity=2).run(testSuites)
    
