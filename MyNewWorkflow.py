import Pipeline as p
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.fasta import Fasta
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.steps.bwa_mem import BwaMem
from Pipeline.bioinformatics.steps.picard_sortsam import PicardSortSam
from Pipeline.bioinformatics.steps.samtools import SamTools

w = p.Workflow("alignment")

inp_bwa_readgroup = p.Input("readgroup", p.String())
inp_reference = p.Input("reference", Fasta())
inp_fastq = p.Input("fastq", Fastq())

inp_gatksortsam_validationstringency = p.Input("validation_stringency", p.String())


step1_bwa = p.Step("bwa", BwaMem())
step2_samtools = p.Step("samtools", SamTools())
step3_gatk_sortsam = p.Step("sortsam", PicardSortSam())

bamOutput = p.Output("outputBam", BamPair())

# all edges for step1
w.add_edges([
    (inp_fastq, step1_bwa.reads),
    (inp_reference, step1_bwa),
    (inp_bwa_readgroup, step1_bwa.readGroup)
])

# step 2
w.add_edge(step1_bwa, step2_samtools)

#step 3
w.add_edges([
    (step2_samtools.out, step3_gatk_sortsam.inputSam),
    (inp_gatksortsam_validationstringency, step3_gatk_sortsam.validation_stringency)
])

w.add_edge(step3_gatk_sortsam.pair, bamOutput)

w.dump_cwl(to_disk=True)