from Pipeline import Workflow, Step, String, Input, Directory, Output
from Pipeline.bioinformatics.data_types.bam import Bam
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline.bioinformatics.tools.bwa.mem.latest import BwaMemLatest
from Pipeline.bioinformatics.tools.gatk4.sortsam.latest import Gatk4SortSamLatest
from Pipeline.bioinformatics.tools.samtools.view.latest import SamToolsViewLatest
from Pipeline.utils.metadata import WorkflowMetadata


class AlignSortedBam(Workflow):

    def __init__(self):
        super(AlignSortedBam, self).__init__("alignsortedbam", friendly_name="Align sorted BAM")

        if not self._metadata:
            self._metadata = WorkflowMetadata()

        self._metadata.documentation = "Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam"
        self._metadata.creator = "Michael Franklin"
        self._metadata.dateCreated = "2018-12-24"
        self._metadata.version = "1.0.0"

        s1_bwa = Step("s1_bwa", BwaMemLatest())
        s2_samtools = Step("s2_samtools", SamToolsViewLatest())
        s3_sortsam = Step("s3_sortsam", Gatk4SortSamLatest())

        s1_inp_header = Input("read_group_header_line", String())
        s1_inp_reference = Input("reference", Fasta())
        s1_inp_fastq = Input("fastq", Fastq())

        s3_inp_tmpdir = Input("tmpdir", Directory())

        o1_bwa = Output("o1_bwa", Sam())
        o2_samtools = Output("o2_samtools", Bam())
        o3_sortsam = Output("o3_sortsam", BamPair())

        # Fully connect step 1
        self.add_edges([
            (s1_inp_header, s1_bwa.readGroupHeaderLine),
            (s1_inp_fastq, s1_bwa.reads),
            (s1_inp_reference, s1_bwa.reference)
        ])
        self.add_default_value(s1_bwa.threads, 36)

        # fully connect step 2
        self.add_edge(s1_bwa, s2_samtools.sam)

        # fully connect step 3
        self.add_edges([
            (s2_samtools.out, s3_sortsam.input),
            (s3_inp_tmpdir, s3_sortsam.tmpDir),
        ])
        self.add_default_value(s3_sortsam.sortOrder, "coordinate")
        self.add_default_value(s3_sortsam.createIndex, True)
        self.add_default_value(s3_sortsam.validationStringency, "SILENT")
        self.add_default_value(s3_sortsam.maxRecordsInRam, 5000000)

        # connect to output
        self.add_edge(s1_bwa, o1_bwa)
        self.add_edge(s2_samtools, o2_samtools)
        self.add_edge(s3_sortsam.output, o3_sortsam)


if __name__ == "__main__":
    print(AlignSortedBam().help())
