import unittest

from Pipeline import Input, String, Step, Directory, Workflow, Array, Output, Int, File
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.bcftools.norm.latest import BcfToolsNormLatest as BcfToolsNorm
from Pipeline.bioinformatics.tools.bwa.mem.latest import BwaMemLatest as BwaMem
from Pipeline.bioinformatics.tools.gatk4.applybqsr.latest import Gatk4ApplyBqsrLatest as Gatk4ApplyBqsr
from Pipeline.bioinformatics.tools.gatk4.baserecalibrator.latest import Gatk4BaseRecalibratorLatest as Gatk4BaseRecalibrator
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.latest import Gatk4GenotypeConcordanceLatest as Gatk4GenotypeConcordance
from Pipeline.bioinformatics.tools.gatk4.haplotypecaller.latest import Gatk4HaplotypeCallerLatest as Gatk4HaplotypeCaller
from Pipeline.bioinformatics.tools.gatk4.markduplicates.latest import Gatk4MarkDuplicatesLatest as Gatk4MarkDuplicates
from Pipeline.bioinformatics.tools.gatk4.mergesamfiles.latest import Gatk4MergeSamFilesLatest as Gatk4MergeSamFiles
from Pipeline.bioinformatics.tools.gatk4.sortsam.latest import Gatk4SortSamLatest as Gatk4SortSam
from Pipeline.bioinformatics.tools.igvtools.index.latest import IgvToolsIndexLatest as IgvToolsIndex
from Pipeline.bioinformatics.tools.samtools.view.latest import SamToolsViewLatest as SamToolsView


class TestGermlinePipeline(unittest.TestCase):


    def create_subworkflow(self):
        sw = Workflow("bwa+st+sort")

        s1_bwa = Step("sw1_bwa", BwaMem())
        s2_samtools = Step("sw_s2", SamToolsView())
        s3_sortsam = Step("sw_s3", Gatk4SortSam())

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", Fasta())
        s1_inp_fastq = Input("fastq", Fastq())

        s3_inp_tmpdir = Input("tmpdir", Directory())

        output = Output("sortedBam", BamPair())

        # Fully connect step 1
        sw.add_edges([
            (s1_inp_header, s1_bwa.readGroupHeaderLine),
            (s1_inp_fastq, s1_bwa.reads),
            (s1_inp_reference, s1_bwa.reference)
        ])
        sw.add_default_value(s1_bwa.threads, 36)

        # fully connect step 2
        sw.add_edge(s1_bwa, s2_samtools.sam)

        # fully connect step 3
        sw.add_edges([
            (s2_samtools.out, s3_sortsam.input),
            (s3_inp_tmpdir, s3_sortsam.tmpDir),
        ])
        sw.add_default_value(s3_sortsam.sortOrder, "coordinate")
        sw.add_default_value(s3_sortsam.createIndex, True)
        sw.add_default_value(s3_sortsam.validationStringency, "SILENT")
        sw.add_default_value(s3_sortsam.maxRecordsInRam, 5000000)

        # connect to output
        sw.add_edge(s3_sortsam.output, output)

        return sw


    def test_workflow(self):

        w = Workflow("whole-genome-germline")

        reference = Input("reference", Fasta())
        fastqInputs = Input("inputs", Array(Fastq()))

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", FastaWithDict())
        s4_inp_SNPS_dbSNP = Input("SNPS_dbSNP", VcfIdx())
        s4_inp_SNPS_1000GP = Input("SNPS_1000GP", VcfIdx())
        s4_inp_OMNI = Input("OMNI", VcfIdx())
        s4_inp_HAPMAP = Input("HAPMAP", VcfIdx())

        inp_tmpdir = Input("tmpdir", Directory())

        s1_sw = Step("bwa+st+sort", self.create_subworkflow())

        s2_mergeSamFiles = Step("s2_mergeSamFiles", Gatk4MergeSamFiles())
        s3_markDuplicates = Step("s3_markDuplicates", Gatk4MarkDuplicates())
        s4_baseRecal = Step("s4_baseRecal", Gatk4BaseRecalibrator())
        s5_applyBqsr = Step("s5_applyBqsr", Gatk4ApplyBqsr())
        s6_haploy = Step("s6_haploy", Gatk4HaplotypeCaller())
        s7_bcfNorm = Step("s7_bcfNorm", BcfToolsNorm())
        s8_bcfNorm = Step("s8_bcfNorm", BcfToolsNorm())
        s9_igv = Step("s9_igv", IgvToolsIndex())
        s10_genotypeConcord = Step("s10_genotypeConcord", Gatk4GenotypeConcordance())


        # step1
        w.add_edges([
            (fastqInputs, s1_sw.fastq),
            (s1_inp_reference, s1_sw.reference),
            (s1_inp_header, s1_sw.read_group_header_line),
            (inp_tmpdir, s1_sw.tmpdir)
        ])

        #step2
        w.add_edge(s1_sw, s2_mergeSamFiles.input)
        w.add_edge(inp_tmpdir, s2_mergeSamFiles.tmpDir)
        w.add_default_value(s2_mergeSamFiles.useThreading, True)
        w.add_default_value(s2_mergeSamFiles.createIndex, True)
        w.add_default_value(s2_mergeSamFiles.maxRecordsInRam, 5000000)
        w.add_default_value(s2_mergeSamFiles.validationStringency, "SILENT")


        # step3
        w.add_edge(s2_mergeSamFiles, s3_markDuplicates)
        w.add_edge(inp_tmpdir, s3_markDuplicates.tmpDir)
        w.add_default_value(s3_markDuplicates.createIndex, True)
        w.add_default_value(s3_markDuplicates.maxRecordsInRam, 5000000)

        #step4 - baserecal
        w.add_edges([
            (s3_markDuplicates, s4_baseRecal),
            (reference, s4_baseRecal.reference),
            (s4_inp_SNPS_dbSNP, s4_baseRecal.knownSites),
            (s4_inp_SNPS_1000GP, s4_baseRecal.knownSites),
            (s4_inp_OMNI, s4_baseRecal.knownSites),
            (s4_inp_HAPMAP, s4_baseRecal.knownSites),
            (inp_tmpdir, s4_baseRecal.tmpDir)
        ])

        # step5 - apply bqsr
        w.add_edges([
            (s3_markDuplicates.output, s5_applyBqsr),
            (s4_baseRecal, s5_applyBqsr),
            (reference, s5_applyBqsr.reference),
            (inp_tmpdir, s5_applyBqsr.tmpDir)
        ])

        # step6 - haplotype caller
        w.add_edges([
            (s5_applyBqsr, s6_haploy),
            (reference, s6_haploy),
            (s4_inp_SNPS_dbSNP, s6_haploy)
        ])

        # step7 - BcfToolsNorm
        # w.add_edges([
        #     ()
        # ])




        o1 = Output("step4_out", File())
        o2_metrics = Output("markDuplicates_metrics", File())

        # connect outputs
        w.add_edges([
            (s3_markDuplicates.metrics, o2_metrics),
            (s4_baseRecal, o1)
        ])

        # w.draw_graph()
        w.dump_cwl(to_disk=True, with_docker=True)


        # temporary_directory = Input("tmpdir", Directory(optional=True))

