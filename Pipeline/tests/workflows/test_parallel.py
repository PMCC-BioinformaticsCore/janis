import unittest

from Pipeline import Workflow, Input, String, Step, Array
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.data_types.vcfidx import VcfIdx
from Pipeline.bioinformatics.steps.bwa_mem import BwaMem
from Pipeline.bioinformatics.steps.gather import Gather
from Pipeline.bioinformatics.steps.gatk_base_recalibrator import GatkRecalibrator
from Pipeline.bioinformatics.steps.gatk_haplotypecaller import GatkHaplotypeCaller
from Pipeline.bioinformatics.steps.gatk_printreads import GatkPrintReads
from Pipeline.bioinformatics.steps.picard_markdup import PicardMarkDup
from Pipeline.bioinformatics.steps.picard_sortsam import PicardSortSam
from Pipeline.bioinformatics.steps.samtools import SamTools


class TestParallel(unittest.TestCase):

    @staticmethod
    def parallel_workflow():
        subworkflow = Workflow("subworkflow")

        sub_inp1_bwa_ref = Input("bwa_ref", Fasta())
        sub_inp1_bwa_reads = Input("reads", Fastq())
        sub_inp1_bwa_read_group = Input("bwa_readGroup", String())

        sub_inp3_sortsam_valstrin = Input("sortsam_validation_stringency", String())

        sub_inp5_recal_ref = Input("recalibrator_ref", FastaWithDict())
        sub_inp5_recal_known = Input("recalibrator_known", Array(VcfIdx()))
        sub_inp5_recal_bed = Input("recalibrator_bed", Bed())



        subworkflow.add_nodes([
            sub_inp1_bwa_read_group,
            sub_inp1_bwa_reads,
            sub_inp1_bwa_ref,

            sub_inp3_sortsam_valstrin,

            sub_inp5_recal_ref,
            sub_inp5_recal_known,
            sub_inp5_recal_bed
        ])

        step1_bwa_mem = Step("bwa_mem", BwaMem())
        step2_samtools= Step("samtools", SamTools())
        step3_sortsam = Step("sortsam", PicardSortSam())
        step4_markdup = Step("markdup", PicardMarkDup())
        # step5_gather = Step("gather", Gather())
        step5_recal = Step("gatk_recal", GatkRecalibrator())
        # step7_printread = Step("gatk_printread", GatkPrintReads())
        # step8_haplo = Step("gatk_haplo", GatkHaplotypeCaller())

        # Add steps
        subworkflow.add_nodes([
            step1_bwa_mem,
            step2_samtools,
            step3_sortsam,
            step4_markdup,
            step5_recal
        ])

        # Inputs -> Step 1
        subworkflow.add_edge(sub_inp1_bwa_ref, step1_bwa_mem)
        subworkflow.add_edge(sub_inp1_bwa_reads, step1_bwa_mem)
        subworkflow.add_edge(sub_inp1_bwa_read_group, step1_bwa_mem)

        # Step1 -> Step 2
        subworkflow.add_edge(step1_bwa_mem, step2_samtools)

        # Step2 -> Step3
        subworkflow.add_edge(step2_samtools, step3_sortsam)
        subworkflow.add_edge(sub_inp3_sortsam_valstrin, step3_sortsam.validation_stringency)

        # Step3 -> Step4
        subworkflow.add_edge(step3_sortsam.out, step4_markdup)

        # Step4 -> Step 5
        subworkflow.add_edge(step4_markdup.outputPair, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_ref, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_known, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_bed, step5_recal)


    def test_parallel(self):


        w = Workflow("parallel")


        sub_normal = self.parallel_workflow()
        sub_tumor = self.parallel_workflow()




