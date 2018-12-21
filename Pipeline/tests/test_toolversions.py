import unittest

from Pipeline.bioinformatics.tools.gatk.haplotypecaller.latest import GatkHaplotypeCallerLatest


class TestToolVersions(unittest.TestCase):


    def test_docker_resolution_order(self):
        tool = GatkHaplotypeCallerLatest()

        # tool.docker() will error if the method order resolution breaks
        self.assertNotEqual(tool.docker(), "")