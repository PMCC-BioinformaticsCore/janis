from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.data_types.sam import Sam
from Pipeline.bioinformatics.data_types.fasta import Fasta

from Pipeline.types.filename import Filename
from Pipeline.types.common_data_types import String, Int, Array, File
from Pipeline.tool.commandtool import CommandTool, ToolOutput, ToolInput, ToolArgument


class BwaMem(CommandTool):

    reads = ToolInput("reads", Fastq(), position=3)
    reference = ToolInput("reference", Fasta(), position=2)
    readGroup = ToolInput("readGroup", String(), prefix="-R")

    minimum_seed_length = ToolInput("minimum_seed_length", Int(optional=True),
                                    position=1, prefix="-k",
                                    doc="-k INT\tminimum seed length")
    threads = ToolInput("threads", Int(optional=True), prefix="-t",
                        default=8, doc="-t INT\tnumber of threads [1]")
    min_std_max_min = ToolInput("min_std_max_min", Array(Int(), optional=True), position=1, prefix="-I")

    outputFilename = ToolInput("outputFilename", Filename(extension='sam'), position=3)

    out = ToolOutput("out", Sam(), glob="$(inputs.outputFilename)")

    @staticmethod
    def tool():
        return "bwa-mem"

    @staticmethod
    def base_command():
        return ["bwa", "mem"]

    @staticmethod
    def docker():
        return "biocontainers/bwa"

    def arguments(self):
        return [ToolArgument("-a", position=2), ToolArgument("-M", position=3)]

    @staticmethod
    def doc():
        return """Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]

  Algorithm options:
         -w INT        band width for banded alignment [100]
         -d INT        off-diagonal X-dropoff [100]
         -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
         -y INT        seed occurrence for the 3rd round seeding [20]
         -c INT        skip seeds with more than INT occurrences [500]
         -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
         -W INT        discard a chain if seeded bases shorter than INT [0]
         -m INT        perform at most INT rounds of mate rescues for each read [50]
         -S            skip mate rescue
         -P            skip pairing; mate rescue performed unless -S also in use
         -e            discard full-length exact matches

  Scoring options:

         -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
         -B INT        penalty for a mismatch [4]
         -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
         -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
         -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
         -U INT        penalty for an unpaired read pair [17]

         -x STR        read type. Setting -x changes multiple parameters unless overriden [null]
                       pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                       ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                       intractg: -B9 -O16 -L5  (intra-species contigs to ref)

  Input/output options:

         -p            smart pairing (ignoring in2.fq)
         -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
         -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
         -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)

         -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
         -T INT        minimum score to output [30]
         -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
         -a            output all alignments for SE or unpaired PE
         -C            append FASTA/FASTQ comment to SAM output
         -V            output the reference FASTA header in the XR tag
         -Y            use soft clipping for supplementary alignments
         -M            mark shorter split hits as secondary

         -I FLOAT[,FLOAT[,INT[,INT]]]
                       specify the mean, standard deviation (10% of the mean if absent), max
                       (4 sigma from the mean if absent) and min of the insert size distribution.
                       FR orientation only. [inferred]

  Note: Please read the man page for detailed description of the command line and options."""


if __name__ == "__main__":
    print(BwaMem().help())