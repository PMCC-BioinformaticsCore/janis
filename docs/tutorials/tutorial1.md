# Tutorial 1 - Building a Workflow

In this stage, we have an installation of Janis, CWLTool, our data and now we're going to construct our workflow. In [Next Generation Sequencing (NGS)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3841808/), short reads of DNA are sequenced in parallel to speed up sequencing time, one of the consequences of many short reads is the need for [alignment](https://en.wikibooks.org/wiki/Next_Generation_Sequencing_(NGS)/Alignment).

After the NGS reads have been pre-processed, we have a `FASTQ` pair, we align the reads inside our FASTQ file into an uncompressed `SAM` file (the _de facto_ standard for short read alignments) using `BWA MEM`, compress this into the binary equivalent `BAM` file using `samtools`, and finally sort the reads using `GATK4 SortSam`.

These tools already exist within the Janis Tool Registry, you can see their documentation right on this website:

- [BWA MEM](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html)
- [Samtols View](https://janis.readthedocs.io/en/latest/tools/bioinformatics/samtools/samtoolsview.html)
- [GATK4 SortSam](https://janis.readthedocs.io/en/latest/tools/bioinformatics/gatk4/gatk4sortsam.html)

## Creating our file

A Janis workflow file is a regular Python file, so we can start by creating a file called `alignment.py` and importing Janis.

```python
import janis
```

## Importing our tools and datatypes

Python requires that you import the tools and types that you use, these import statements are available on the documentation. We'll have one import per tool, and one import for every data-types we use.

We have three inputs we want to expose on this workflow:

1. Sample name (`String`)
2. Sequencing Reads (`FastqGzPair` - paired end sequence)
3. Reference files (`Fasta` + index files)


We can use the `janis.String` datatype (imported with Janis) for the first, and we can import the remaining bioinformatics types:

```python
from janis_bioinformatics.tools.bwa.mem.latest import BwaMemLatest
from janis_bioinformatics.tools.samtools.view.view import SamToolsView_1_9
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_1_2
from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict
```

## Declaring our workflow and exposing inputs

We'll create an instance of the [`janis.WorkflowBuilder`](https://janis.readthedocs.io/en/latest/references/workflow.html#janis.Workflow) class, this requires a workflow identifier. We discussed which imports we want in the previous section which we can expose on this workflow with the `janis.Workflow.input` method:

```python
w = janis.WorkflowBuilder("alignmentWorkflow")

w.input("readGroup", janis.String)
w.input("fastq", FastqGzPair)
w.input("reference", FastaWithDict)
```

## Declaring our steps and connections

Steps are easy to create, however you may need to refer to the documentation when writing your own workflows to know which named parameters a tool takes (and their types). Similar to exposing inputs, we create steps with the `janis.Workflow.step` method.

We can refer to any node on the workflow graph (such as an input) by accessing the property of the same name (dot-notation). Eg, to access the `sampleName` on our workflow, we can use `w.sampleName`.

We instantiate our tool with the named parameters we want to provide and pass that as the second parameter to the `janis.Workflow.step` method:

### BWA MEM

```python
w.step(
    "bwamem", 
    BwaMemLatest( 
        reads=w.fastq, 
        sampleName=w.sampleName, 
        reference=w.reference
    )
)
```

### Samtools view

When creating the connection between `bwamem` and `samtoolsview`, we'll access the `out` output of `BwaMemLatest`. This will create a dependency of `"bwamem"` for `samtoolsview`.

```python
w.step("samtoolsview", SamToolsView_1_9(sam=w.bwamem.out))
```

### SortSam

SortSam requires a number of values we want to set 

```python
w.step(
    "sortsam",
    Gatk4SortSam_4_1_2(
        bam=w.samtoolsview.out,
        sortOrder="coordinate",
        createIndex=True,
        validationStringency="SILENT",
        maxRecordsInRam=5000000
    )
)
```

## Exposing outputs

Outputs have a very similar syntax to both inputs and steps, they take an `identifier` and a named `source` parameter. We only want to output the resulting bam file from `sortsam`, which we can do with the following line:

```python
w.output("out", source=w.sortsam.out)
```

## Workflow + Translation

Hopefully you have a workflow that looks like the following!

```python
import janis

from janis_bioinformatics.tools.bwa.mem.latest import BwaMemLatest
from janis_bioinformatics.tools.samtools.view.view import SamToolsView_1_9
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_1_2
from janis_bioinformatics.data_types import FastqGzPair, FastaWithDict

w = janis.WorkflowBuilder("alignmentWorkflow")

# Inputs
w.input("sampleName", janis.String, value="sampleName")
w.input("fastq", FastqGzPair, value="/path/to/reads.fastq")
w.input("reference", FastaWithDict, value="/path/to/reference.fasta")

# Steps
w.step(
    "bwamem", 
    BwaMemLatest( 
        reads=w.fastq, 
        sampleName=w.sampleName, 
        reference=w.reference
    )
)
w.step("samtoolsview", SamToolsView_1_9(sam=w.bwamem.out))
w.step(
    "sortsam",
    Gatk4SortSam_4_1_2(
        bam=w.samtoolsview.out,
        sortOrder="coordinate",
        createIndex=True,
        validationStringency="SILENT",
        maxRecordsInRam=5000000
    )
)

# Outputs
w.output("out", source=w.sortsam.out)
```

We can translate the following file into Workflow Description Language using janis from the terminal:

```bash
janis translate alignment.py wdl
```
