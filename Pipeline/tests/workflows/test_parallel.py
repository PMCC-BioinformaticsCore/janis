import unittest

from Pipeline import Workflow, Input, String, Step
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.steps.bwa_mem import BwaMem
from Pipeline.bioinformatics.steps.picard_sortsam import PicardSortSam
from Pipeline.bioinformatics.steps.samtools import SamTools


class TestParallel(unittest.TestCase):

    @staticmethod
    def parallel_workflow():
        subworkflow = Workflow("subworkflow")

        sub_inp_bwa_ref = Input("bwa_ref", Fasta())
        sub_inp_bwa_reads = Input("reads", Fastq())
        sub_inp_bwa_read_group = Input("bwa_readGroup", String())

        sub_inp_sortsam_valstrin = Input("sortsam_validation_stringency", String())

        step1_bwa_mem = Step("bwa_mem", BwaMem())
        step2_samtools= Step("samtools", SamTools())
        step3_sortsam = Step("sortsam", PicardSortSam())

        subworkflow.add_edge(sub_inp_bwa_ref, step1_bwa_mem)
        subworkflow.add_edge(sub_inp_bwa_reads, step1_bwa_mem)

        subworkflow.add_edge(sub_inp_bwa_read_group, step2_samtools)
        subworkflow.add_edge(step1_bwa_mem, step2_samtools)

        subworkflow.add_edge(sub_inp_sortsam_valstrin, step3_sortsam)
        subworkflow.add_edge(step2_samtools, step3_sortsam)


    def test_parallel(self):


        w = Workflow("parallel")


        sub_normal = self.parallel_workflow()
        sub_tumor = self.parallel_workflow()




