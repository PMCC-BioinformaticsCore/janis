# Tutorial 3 - Naming and organising outputs

Sometimes it's useful to organise your outputs in specific ways, especially when Janis might need to interact with other specific systems. By default, outputs are named by the tag of the output. The extension is derived from the type of the output (the [Bam](https://janis.readthedocs.io/en/latest/datatypes/bam.html) filetype knows it uses a `.bam` extension).

For example, if you had an output called `out` which had type `BamBai`, your output files by default would have the name `out.bam` and `out.bam.bai`.

When exposing an output on a workflow, there are two arguments you can provide to override this behaviour:

- `output_name: Union[str, InputSelector, InputNode] = None`
- `output_folder: Union[str, Tuple[Node, str], List[Union[str, InputSelector, Tuple[Node, str]]]] = None`

## Preparation

This tutorial uses the workflow build in [Tutorial 1](https://janis.readthedocs.io/en/latest/tutorials/tutorial1.html). You can follow this guide, but  the example will require you to have built (or acquired) the workflow from tutorial 1.


## Output name

Simply put, `output_name` is the dervied filename of the output without the extension. By default, this is the `tag` of the output.

You can specify a new output name in 2 ways:

1. A static string: `output_name="new name for output"`
2. Selecting an input value (given "sample_name" is the name of an input):
    1. `output_name=workflow.sample_name`
    2. `output_name=InputSelector("sample_name")`

You should make the following considerations:

- The input you select should be a string, or

- If the output you're naming is an array, the input you select should either be:
    - singular
    - have the same number of elements in it.

    Janis will either fall back to the first element if it's a list, or default to the output tag. This may cause outputs to override each other.


## Output folder

Similar to the output name, the `output_folder` is folder, or group of nested folders into which your output will be written. By default, this field has no value and outputs are linked directly into the output directory.

If the output_folder field is an array, a nested folder is created for each element in ascending order (eg: `["parent", "child", "child_of_child"]`).

There are multiple ways to specify output directories:

1. A static string: `output_folder="my_outputs"`
2. Selecting an input value (given "sample_name" is the name of an input):
    1. `output_folder=workflow.sample_name`
    2. `output_folder=InputSelector("sample_name")`
3. An array of a combination of values:
    - `output_folder=["variants", "unmerged"]`
    - `output_folder=["variants", w.sample_name]`
    - `output_folder=[w.other_field, w.sample_name]`

## Alignment workflow

In our alignment example (`tools/alignment.py`), we have the following outputs:

```python
w.output("out", source=w.sortsam.out)
w.output("stats", source=w.flagstat.stats)
```

We want to name the output in the following way:

- Grouped the bam `out` into the folder: `bams`,
- Group the samtools summary into folder `statistics`
- Both outputs named: `sample_name` (Janis will automatically add the `.bam` or `.txt` extension).

```python
w.output(
    "out", 
    source=w.sortsam.out,
    output_name=w.sample_name,
    output_folder=["bams"]
)
w.output(
    "stats", 
    source=w.flagstat.stats,
    output_name=w.sample_name,
    output_folder=["statistics"]
)
```

Run the workflow again and inspect the new results:

```bash
janis translate tools/alignment.py wdl
janis run -o tutorial3 tools/alignment.py \
    --fastq data/BRCA1_R*.fastq.gz \
    --reference reference/hg38-brca1.fasta \
    --sample_name NA12878 \
    --read_group "@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPL:ILLUMINA"
```

And then after the workflow has finished, we find the following output structure:

```
$ ls tutorial3/*

tutorial3/bams:
-rw-r--r--  PMCI\Bioinf-Cluster   2.7M NA12878.bam
-rw-r--r--  PMCI\Bioinf-Cluster   296B NA12878.bam.bai

tutorial3/statistics:
-rw-r--r--  PMCI\Bioinf-Cluster   408B NA12878.txt

tutorial3/janis:
[...ommitted]
```


