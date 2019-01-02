import unittest

from Pipeline import Input, String, Step, Directory, Workflow, Array, Output, Int
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastq import Fastq
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

        s1 = Step("s1", BwaMem())
        s2 = Step("sw_s2", SamToolsView())
        s3 = Step("sw_s3", Gatk4SortSam())

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", Fasta())
        s1_inp_fastq = Input("fastq", Fastq())

        s3_inp_tmpdir = Input("tmpdir", Directory())
        s3_inp_sortorder = Input("sortorder", String()) #default="coordinate"))
        s3_inp_maxramrecords = Input("max_records_in_ram", Int())#default=5000000))
        s3_validation_stringency = Input("validation_stringency", String())#default="SILENT"))

        output = Output("sortedBam", Bam())


        # Fully connect step 1
        sw.add_edges([
            (s1_inp_header, s1.readGroupHeaderLine),
            (s1_inp_fastq, s1.reads),
            (s1_inp_reference, s1.reference)
        ])

        # fully connect step 2
        sw.add_edges([
            (s1, s2.sam)
        ])

        # fully connect step 3
        sw.add_edges([
            (s2.out, s3.input),
            (s3_inp_tmpdir, s3.tmpDir),
            (s3_inp_sortorder, s3.sortOrder),
            (s3_inp_maxramrecords, s3.maxRecordsInRam),
            (s3_validation_stringency, s3.validationStringency)
        ])

        # connect to output
        sw.add_edge(s3.output, output)

        return sw


    def test_workflow(self):

        w = Workflow("whole-genome-germline")

        reference = Input("reference", Fasta())
        fastqInputs = Input("inputs", Array(Fastq()))

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", Fasta())

        s1_inp_tmpdir = Input("tmpdir", Directory())
        s1_inp_sortorder = Input("sortorder", String()) #default="coordinate"))
        s1_inp_maxramrecords = Input("max_records_in_ram", Int())#default=5000000))
        s1_validation_stringency = Input("validation_stringency", String())#default="SILENT"))

        s1_sw = Step("bwa+st+sort", self.create_subworkflow())

        s2 = Step("s2", Gatk4MergeSamFiles())
        s3 = Step("s3", Gatk4MarkDuplicates())
        s4 = Step("s4", Gatk4BaseRecalibrator())
        s5 = Step("s5", Gatk4ApplyBqsr())
        s6 = Step("s6", Gatk4HaplotypeCaller())
        s7 = Step("s7", BcfToolsNorm())
        s8 = Step("s8", BcfToolsNorm())
        s9 = Step("s9", IgvToolsIndex())
        s10 = Step("s10", Gatk4GenotypeConcordance())


        # step1
        w.add_edges([
            (fastqInputs, s1_sw.fastq),
            (s1_inp_reference, s1_sw.reference),
            (s1_inp_header, s1_sw.read_group_header_line),
            (s1_inp_tmpdir, s1_sw.tmpdir),
            (s1_inp_sortorder, s1_sw.sortorder),
            (s1_inp_maxramrecords, s1_sw.max_records_in_ram),
            (s1_validation_stringency, s1_sw.validation_stringency)
        ])

        #step2
        w.add_edges([
            (s1_sw, s2)
        ])

        #step3
        w.add_edges([
            (s2, s3)
        ])

        #step4
        w.add_edges([
            (s3, s4),
            (reference, s3.reference),
            (snps_dbsnnp, step4.knownSites)
        ])

        w.add_edges([
            (s4, s5)
        ])

        o1 = Output( "step4_out", File())

        w.add_edge(s4, o1)


        temporary_directory = Input("tmpdir", Directory(optional=True))

