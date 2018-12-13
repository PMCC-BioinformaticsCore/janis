import unittest

from Pipeline import Workflow, Input, String, Step
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.steps.bwa_mem import BwaMem
from Pipeline.bioinformatics.steps.samtools import SamTools


class TestParallel(unittest.TestCase):

    @staticmethod
    def parallel_workflow():
        subworkflow = Workflow("subworkflow")

        sub_inp_bwa_ref = Input("bwa_ref", Fasta())
        sub_inp_bwa_reads = Input("reads", Fastq())
        sub_inp_bwa_read_group = Input("bwa_readGroup", String())

        step1_bwa_mem = Step("bwa_mem", BwaMem())

        step2_samtools= Step("samtools", SamTools())



    def test_parallel(self):


        w = Workflow("parallel")


        sub_normal = self.parallel_workflow()
        sub_tumor = self.parallel_workflow()




