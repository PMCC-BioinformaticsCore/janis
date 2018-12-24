from Pipeline import Workflow, Input, Step
from Pipeline.bioinformatics.data_types.fastq import Fastq
from Pipeline.bioinformatics.steps.bwa_mem import BwaMem
from Pipeline.bioinformatics.steps.fastqc import FastQC


def test_workflow():
    w = Workflow("amit")

    fastq = Input("fastq", Fastq())

    step1 = Step("fastqc", FastQC())
    step2 = Step("bwa", BwaMem())


    w.add_edges([
        (fastq, step1),
        (step1, step2)
    ])

    # w.add_pipe(fastq, step1, step2)

    w.dump_cwl(to_disk=False)




