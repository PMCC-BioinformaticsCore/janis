import unittest

from Pipeline import Input, String, Step, Directory, Workflow, Array, Output, Int, File
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastawithdict import FastaWithDict
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline.bioinformatics.data_types.vcf import VcfIdx, TabixIdx
from Pipeline.bioinformatics.tools.bcftools.norm.latest import BcfToolsNormLatest as BcfToolsNorm
from Pipeline.bioinformatics.tools.bwa.mem.latest import BwaMemLatest as BwaMem
from Pipeline.bioinformatics.tools.gatk4.applybqsr.latest import Gatk4ApplyBqsrLatest as Gatk4ApplyBqsr
from Pipeline.bioinformatics.tools.gatk4.baserecalibrator.latest import Gatk4BaseRecalibratorLatest as Gatk4BaseRecalibrator
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.latest import Gatk4GenotypeConcordanceLatest as Gatk4GenotypeConcordance
from Pipeline.bioinformatics.tools.gatk4.haplotypecaller.latest import Gatk4HaplotypeCallerLatest as Gatk4HaplotypeCaller
from Pipeline.bioinformatics.tools.gatk4.markduplicates.latest import Gatk4MarkDuplicatesLatest as Gatk4MarkDuplicates
from Pipeline.bioinformatics.tools.gatk4.mergesamfiles.latest import Gatk4MergeSamFilesLatest as Gatk4MergeSamFiles
from Pipeline.bioinformatics.tools.gatk4.sortsam.latest import Gatk4SortSamLatest as Gatk4SortSam
from Pipeline.bioinformatics.tools.htslib.bgzip.latest import BGZipLatest
from Pipeline.bioinformatics.tools.htslib.tabix.latest import TabixLatest
from Pipeline.bioinformatics.tools.igvtools.index.latest import IgvToolsIndexLatest as IgvToolsIndex
from Pipeline.bioinformatics.tools.samtools.view.latest import SamToolsViewLatest as SamToolsView


class TestGermlinePipeline(unittest.TestCase):


    def create_subworkflow(self):
        sw = Workflow("bwa_st_sort")

        s1_bwa = Step("s1_bwa", BwaMem())
        s2_samtools = Step("s2_samtools", SamToolsView())
        s3_sortsam = Step("s3_sortsam", Gatk4SortSam())

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", Fasta())
        s1_inp_fastq = Input("fastq", Fastq())

        s3_inp_tmpdir = Input("tmpdir", Directory())

        o1_bwa = Output("o1_bwa", Sam())
        o2_samtools = Output("o2_samtools", Bam())
        o3_sortsam = Output("o3_sortsam", BamPair())

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
        sw.add_edge(s1_bwa, o1_bwa)
        sw.add_edge(s2_samtools, o2_samtools)
        sw.add_edge(s3_sortsam.output, o3_sortsam)

        return sw

    def test_workflow(self):

        w = Workflow("whole-genome-germline")

        reference = Input("reference", Fasta())
        fastqInputs = Input("inputs", Array(Fastq()))

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", FastaWithDict())
        s4_inp_SNPS_dbSNP = Input("SNPS_dbSNP", VcfIdx())
        s4_inp_SNPS_1000GP = Input("SNPS_1000GP", TabixIdx())
        s4_inp_OMNI = Input("OMNI", TabixIdx())
        s4_inp_HAPMAP = Input("HAPMAP", TabixIdx())
        s10_truth = Input("TRUTH_VCF", VcfIdx())
        s10_intervals = Input("INTERVALS", Array(VcfIdx()))

        inp_tmpdir = Input("tmpdir", Directory())

        s1_sw = Step("sw_bwa_st_sort", self.create_subworkflow())

        s2_mergeSamFiles = Step("s2_mergeSamFiles", Gatk4MergeSamFiles())
        s3_markDuplicates = Step("s3_markDuplicates", Gatk4MarkDuplicates())
        s4_baseRecal = Step("s4_baseRecal", Gatk4BaseRecalibrator())
        s5_applyBqsr = Step("s5_applyBqsr", Gatk4ApplyBqsr())
        s6_haploy = Step("s6_haploy", Gatk4HaplotypeCaller())
        s7_bcfNorm = Step("s7_bcfNorm", BcfToolsNorm())
        s8_bgzip = Step("s8_bgzip", BGZipLatest())
        s9_tabix = Step("s9_tabix", TabixLatest())
        s10_genotypeConcord = Step("s10_genotypeConcord", Gatk4GenotypeConcordance())

        # step1
        w.add_edge(fastqInputs, s1_sw.fastq)
        w.add_edges([
            (s1_inp_reference, s1_sw.reference),
            (s1_inp_header, s1_sw.read_group_header_line),
            (inp_tmpdir, s1_sw.tmpdir)
        ])

        #step2
        w.add_edge(s1_sw.o3_sortsam, s2_mergeSamFiles.input)
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
        w.add_edges([
            (reference, s7_bcfNorm),
            (s6_haploy, s7_bcfNorm)
        ])

        # step8 - BGZip
        w.add_edge(s7_bcfNorm, s8_bgzip.file)

        # step9 - tabix
        w.add_edge(s8_bgzip, s9_tabix.file)

        # step10 - gatk-genotypeconcordance
        w.add_edges([
            (s9_tabix, s10_genotypeConcord.callVCF),
            (s10_truth, s10_genotypeConcord.truthVCF),
            (s10_intervals, s10_genotypeConcord.intervals),
        ])
        w.add_default_value(s10_genotypeConcord.treatMissingSitesAsHomeRef, True)


        # Outputs

        w.add_edges([
            (s1_sw.o1_bwa, Output("sw_bwa")),
            (s1_sw.o2_samtools, Output("sw_samtools")),
            (s1_sw.o3_sortsam, Output("sw_sortsam")),
            (s2_mergeSamFiles, Output("o4_merged")),
            (s3_markDuplicates.output, Output("o5_marked_output")),
            (s3_markDuplicates.metrics, Output("o5_marked_metrics")),
            (s4_baseRecal, Output("o6_recal")),
            (s5_applyBqsr, Output("o7_bqsr")),
            (s6_haploy.output, Output("o8_halpo")),
            (s7_bcfNorm, Output("o9_bcfnorm")),
            (s8_bgzip, Output("o10_bgzip")),
            (s9_tabix, Output("o11_tabix")),
            (s10_genotypeConcord.vcf, Output("o12_concord")),
            (s10_genotypeConcord.summaryMetrics, Output("o12_concord")),
            (s10_genotypeConcord.detailMetrics, Output("o12_concord")),
            (s10_genotypeConcord.contingencyMetrics, Output("o12_concord"))
        ])

        # w.draw_graph()
        w.dump_cwl(to_disk=True, with_docker=False)


        # temporary_directory = Input("tmpdir", Directory(optional=True))

