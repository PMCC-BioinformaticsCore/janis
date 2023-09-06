# Tutorial 1 - Building a Workflow

In this stage, we're going to build a simple workflow to align short reads of DNA.

1. Start with a pair of compressed `FASTQ` files,
2. Align these reads using `BWA MEM` into an uncompressed `SAM` file (the _de facto_ standard for short read alignments),
3. Compress this into the binary equivalent `BAM` file using `samtools`, and finally
4. Sort the reads using `GATK4 SortSam`.

These tools already exist within the Janis Tool Registry, you can see their documentation online:

- [BWA MEM](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html)
- [Samtols View](https://janis.readthedocs.io/en/latest/tools/bioinformatics/samtools/samtoolsview.html)
- [GATK4 SortSam](https://janis.readthedocs.io/en/latest/tools/bioinformatics/gatk4/gatk4sortsam.html)

## Preparation

To prepare for this tutorial, we're going to create a folder and download some data:

```bash
mkdir janis-tutorials && cd janis-tutorials

# If WGET is installed
wget -q -O- "https://github.com/PMCC-BioinformaticsCore/janis-workshops/raw/master/janis-data.tar" | tar -xz

# If CURL is installed
curl -Ls "https://github.com/PMCC-BioinformaticsCore/janis-workshops/raw/master/janis-data.tar" | tar -xz
```


## Creating our file

A Janis workflow is a Python script, so we can start by creating a file called `alignment.py` and importing Janis.

```bash
mkdir tools
vim tools/alignment.py # or vim, emacs, sublime, vscode
```

From the `janis_core` library, we're going to import `WorkflowBuilder` and a `String`:

```python
from janis_core import WorkflowBuilder, String
```

## Imports

We have four inputs we want to expose on this workflow:

1. Sequencing Reads (`FastqGzPair` - paired end sequence)
2. Sample name (`String`)
3. Read group header (`String`)
4. Reference files (`Fasta` + index files (`FastaWithIndex`))

We've already imported the `String` type, and we can import `FastqGzPair` and `FastaWithIndex` from the `janis_bioinformatics` registry:

```python
from janis_bioinformatics.data_types import FastqGzPairedEnd, FastaWithDict
```

### Tools

We've discussed the tools we're going to use. The documentation for each tool has a row in the tbale caled "Python" that gives you the import statement. This is how we'll import how tools:


```python
from janis_bioinformatics.tools.bwa import BwaMemLatest
from janis_bioinformatics.tools.samtools import SamToolsView_1_9
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_1_2
```



## Declaring our workflow

We'll create an instance of the [`WorkflowBuilder`](https://janis.readthedocs.io/en/latest/references/workflow.html#janis.Workflow) class, this just requires a name for your workflow (can contain alphanumeric characters and underscores).

```python
w = WorkflowBuilder("alignmentWorkflow")
```

A workflow has 3 methods for building workflows:

- `workflow.input` - Used for creating inputs,
- `workflow.step` - Creates a step on a workflow,
- `workflow.output` - Exposes an output on a workflow.

We give each input / step / output a unique identifier, which then becomes a node in our workflow graph. We can refer to the created node using _dot-notation_ (eg: `w.input_name`). We'll see how this works in the later sections.

More information about each step will be linked from this page about the [`Workflow` and `WorkflowBuilder` class](https://janis.readthedocs.io/en/latest/references/workflow.html).


### Creating inputs on a workflow

> Further reading: [Creating an input](https://janis.readthedocs.io/en/latest/references/workflow.html#creating-an-input)

To create an input on a workflow, you can use the `Workflow.input` method, which has the following structure:

```python
Workflow.input(
    identifier: str, 
    datatype: DataType, 
    default: any = None, 
    doc: str = None
)
```

An input requires a unique identifier (string) and a DataType (String, FastqGzPair, etc). Let's prepare the inputs for our workflow:


```python
w.input("sample_name", String)
w.input("read_group", String)
w.input("fastq", FastqGzPairedEnd)
w.input("reference", FastaWithDict)
```

### Declaring our steps and connections

> Further reading: [Creating a step](https://janis.readthedocs.io/en/latest/references/workflow.html#creating-a-step)

Similar to exposing inputs, we create steps with the `Workflow.step` method. It has the following structure:

```python
Workflow.step(
    identifier: str, 
    tool: janis_core.tool.tool.Tool, 
    scatter: Union[str, List[str], ScatterDescription] = None, 
)
```

We provide a identifier for the step (unique amongst the other nodes in the workflow), and intialise our tool, passing our inputs of the step as parameters.

We can refer to an input (or previous result) using the dot notation. For example, to refer to the `fastq` input, we can use `w.fastq`.

#### BWA MEM

We use [bwa mem's documentation](https://janis.readthedocs.io/en/latest/tools/bioinformatics/bwa/bwamem.html) to determine that we need to provide the following inputs:

- `reads`: `FastqGzPair`            (connect to `w.fastq`)
- `readGroupHeaderLine`: `String`   (connect to `w.read_group`)
- `reference`: `FastaWithDict`      (connect to `w.reference`)

We can connect them to the relevant inputs to get the following step definition:

```python
w.step(
    "bwamem",   # identifier
    BwaMemLatest(
        reads=w.fastq,
        readGroupHeaderLine=w.read_group,
        reference=w.reference
    )
)
```

#### Samtools view

We'll use a very similar pattern for Samtools View, except this time we'll reference the output of `bwamem`. From bwa mem's documentation, there is one output called `out` with type `Sam`. We'll connect this to `SamtoolsView` only input, called `sam`.


```python
w.step(
    "samtoolsview",
    SamToolsView_1_9(
        sam=w.bwamem.out
    )
)
```

#### SortSam

In addition to connecting the output of `samtoolsview` to Gatk4 SortSam, we want to tell SortSam to use the following values:

- sortOrder: `"coordinate"`
- createIndex: `True`

Instead of connecting an input or a step, we just just provide the literal value.

```python
w.step(
    "sortsam",
    Gatk4SortSam_4_1_2(
        bam=w.samtoolsview.out,
        sortOrder="coordinate",
        createIndex=True
    )
)
```

### Exposing outputs

> Further reading: [Creating an output](https://janis.readthedocs.io/en/latest/references/workflow.html#creating-an-output)

Outputs have a very similar syntax to both inputs and steps, they take an `identifier` and a named `source` parameter. Here is the structure:

```python
Workflow.output(
    identifier: str,
    datatype: DataType = None,
    source: Node = None,
    output_folder: List[Union[String, Node]] = None,
    output_name: Union[String, Node] = None
)
```

Often, we don't want to specify the output data type, because we can let Janis do this for us. We'll talk about the `output_folder` and `output_name` in the next few sections. For now, we just have to specify an output identifier and a source.

```python
w.output("out", source=w.sortsam.out)
```

## Workflow + Translation

Hopefully you have a workflow that looks like the following!

```python
from janis_core import WorkflowBuilder, String

from janis_bioinformatics.data_types import FastqGzPairedEnd, FastaWithDict

from janis_bioinformatics.tools.bwa import BwaMemLatest
from janis_bioinformatics.tools.samtools import SamToolsView_1_9
from janis_bioinformatics.tools.gatk4 import Gatk4SortSam_4_1_2

w = WorkflowBuilder("alignmentWorkflow")

# Inputs
w.input("sample_name", String)
w.input("read_group", String)
w.input("fastq", FastqGzPairedEnd)
w.input("reference", FastaWithDict)

# Steps
w.step(
    "bwamem", 
    BwaMemLatest( 
        reads=w.fastq, 
        readGroupHeaderLine=w.read_group, 
        reference=w.reference
    )
)
w.step(
    "samtoolsview", 
    SamToolsView_1_9(
        sam=w.bwamem.out
    )
)

w.step(
    "sortsam",
    Gatk4SortSam_4_1_2(
        bam=w.samtoolsview.out,
        sortOrder="coordinate",
        createIndex=True
    )
)

# Outputs
w.output("out", source=w.sortsam.out)
```

We can translate the following file into Workflow Description Language using janis from the terminal:

```bash
janis translate tools/alignment.py wdl
```


## Running the alignment workflow

We'll run the workflow against the current directory.

```bash
janis run -o . --engine cwltool \
    tools/alignment.py \
    --fastq data/BRCA1_R*.fastq.gz \
    --reference reference/hg38-brca1.fasta \
    --sample_name NA12878 \
    --read_group "@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPL:ILLUMINA"
```

After the workflow has run, you'll see the outputs in the current directory:

```bash
ls

# drwxr-xr-x  mfranklin  1677682026   160B  data
# drwxr-xr-x  mfranklin  1677682026   256B  janis
# -rw-r--r--  mfranklin  wheel        2.7M  out.bam
# -rw-r--r--  mfranklin  wheel        296B  out.bam.bai
# drwxr-xr-x  mfranklin  1677682026   320B  reference
# drwxr-xr-x  mfranklin  1677682026   128B  tools
```

### OPTIONAL: Run with Cromwell

If you have `java` installed, Janis can run the workflow in the Crowmell execution engine by using the `--engine cromwell` parameter:

```bash
janis run -o run-with-cromwell --engine cromwell \
    tools/alignment.py \
    --fastq data/BRCA1_R*.fastq.gz \
    --reference reference/hg38-brca1.fasta \
    --sample_name NA12878 \
    --read_group "@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPL:ILLUMINA"
```
