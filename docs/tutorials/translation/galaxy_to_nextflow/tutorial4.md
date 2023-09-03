
# Tutorial 4: Simple Galaxy Tool -> Nextflow

<br>

## Sections

- [Introduction](#introduction)
- [Running Janis Translate](#running-janis-translate)
- [Manual Adjustments](#manual-adjustments)
- [Running Samtools Flagstat as a Workflow](#running-samtools-flagstat-as-a-workflow)
- [Conclusion](#conclusion)

<br>

## Introduction

This section demonstrates translation of a basic `samtools flagstat` Galaxy Tool Wrapper to Nextflow using `janis translate`. <br>

The Galaxy Tool Wrapper used in this section was created by contributors to the [Galaxy Devteam repository](https://github.com/galaxyproject/tools-devteam) of tools. 

The underlying software used in the Galaxy Tool Wrapper - [samtools_flagstat](http://www.htslib.org/doc/samtools-flagstat.html) - displays summary information for an alignment file. 

<br>

**Software**

Before continuing, ensure you have the following software installed:
- [Nextflow](https://nf-co.re/usage/installation)
- [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) or [Docker](https://docs.docker.com/engine/install/)
- [Janis](https://janis.readthedocs.io/en/latest/index.html)

<br>

**IDE**

Any IDE or a CLI text editor (VIM, nano) are sufficient for this material. 

We recommend Visual Studio Code (VS Code) as it is lightweight and has rich support for extensions to add functionality. 

<br>

**Obtaining Janis**

In this tutorial we will use a singularity container to run `janis translate`. 

Containers are great because they remove the need for package managers, and guarantee that the software can run on any machine. 

Run the following command to pull the janis image:
```
singularity pull janis.sif docker://pppjanistranslate/janis-translate:0.13.0
```

Check your image by running the following command:
```
singularity exec ~/janis.sif janis translate
```

If the image is working, you should see the janis translate helptext.

<br>

**Downloading Source Files and Sample Data**

For this tutorial we will fetch all necessary data from zenodo using wget.  

This archive contains sample data and the finished translations as a reference.

Run the following commands to download & decompress the zenodo archive:
```
wget https://zenodo.org/record/8275567/files/tutorial4.tar.gz
tar -xvf tutorial4.tar.gz
```

After you have decompressed the tar archive, change into the new directory: 

```
cd tutorial4
```

Inside this folder we have the following structure: 

```
tutorial4
├── data
│   └── alignments.bam
└── final
    ├── nextflow.config
    └── samtools_flagstat.nf
``` 

We will translate the Galaxy "samtools flagstat" Tool Wrapper into nextflow using janis, then will test out our translation using the bam file in the `data/` folder. 

<br>

## Running Janis Translate

To translate a tool / workflow,  we use `janis translate`.

```
janis translate --from <src> --to <dest> <filepath>
```

The `--from` argument specifies the workflow language of the source file(s), and `--to` specifies the destination we want to translate to. 

In our case, this will be `--from galaxy --to nextflow`.

The `<filepath>` argument is the source file we will translate. 

Aside from local filepaths, `janis translate` can also access Galaxy Tool Wrappers using a tool ID. 
We will use this method today as it is an easier way to access Tool Wrappers. 

To get the `samtools_flagstat` Tool ID, navigate to the tool using any usegalaxy.org server. <br>
The following is a link to the samtools flagstat tool in Galaxy Australia: <br> 
https://usegalaxy.org.au/root?tool_id=toolshed.g2.bx.psu.edu/repos/devteam/samtools_flagstat/samtools_flagstat/2.0.4. 

Once here, we will copy the Tool ID. 

![alt text](../media/samtools_flagstat_copy_tool_id.PNG)

At time of writing, the current Tool ID for the `samtools_flagstat` tool wrapper is *toolshed.g2.bx.psu.edu/repos/devteam/samtools_flagstat/samtools_flagstat/2.0.4*

Now we have the Tool ID, we can access & translate this Galaxy Tool Wrapper to a Nextflow process. 

To translate the Galaxy Tool, run the following command:
```
singularity exec ~/janis.sif janis translate -o samtools_flagstat --from galaxy --to nextflow toolshed.g2.bx.psu.edu/repos/devteam/samtools_flagstat/samtools_flagstat/2.0.4
```

<br>

Once complete, you will see a folder called `translated/` appear, and a nextflow process called `samtools_flagstat.nf` will be present inside. 

For your own reference / interest, the actual Galaxy Tool Wrapper files will be downloaded during translation & will be presented to you in `translated/source/`. 

<br>

## Manual Adjustments

The `translated/samtools_flagstat.nf` file should be similar to the following: 

```
nextflow.enable.dsl=2

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.13--h8c37831_0"

    input:
    path input1
    val addthreads

    output:
    path "output1.txt", emit: output1

    script:
    """
    samtools flagstat \
    --output-fmt "txt" \
    -@ ${addthreads} \
    ${input1} \
    > output1.txt
    """

}
```

We can see that this nextflow process has two inputs, a single output, and calls `samtools flagstat`.  

Before continuing, let's check the [samtools flagstat](http://www.htslib.org/doc/samtools-flagstat.html) documentation. 
In the documentation, we see the following:
```
samtools flagstat in.sam|in.bam|in.cram

-@ INT
Set number of additional threads to use when reading the file.

--output-fmt/-O FORMAT
Set the output format. FORMAT can be set to `default', `json' or `tsv' to select the default, JSON or tab-separated values output format. If this option is not used, the default format will be selected.
```

By matching up the process `inputs:` section and the `script:` section, we can see that:
- `path input1` will be the input `sam | bam | cram`
- `val addthreads` will be the threads argument passed to `-@`
- the `--output-fmt` option has been assigned the default value of `"txt"`

We can also see that a container image is available for this tool. 

This translation is correct for the `samtools flagstat` Galaxy tool wrapper and needs no adjusting. 

<br>

> Note: <br>
> If you would like to expose the `--output-fmt` option as a process input, you can do the following: 
>
> - add a `val format` input to the process 
> - reference this input in the script, replacing the hardcoded `"txt"` value<br>
>   (e.g. `--output-fmt ${format}`)

<br>

## Running Samtools Flagstat as a Workflow

**Setting up nextflow.config**

To run this process, we will set up a `nextflow.config` file and add some lines to the top of our process definition to turn it into a workflow.

Create a new file called `nextflow.config` in the `translated/` folder alongside `samtools_flagstat.nf`. 

Copy and paste the following code into your `nextflow.config` file: 

```
nextflow.enable.dsl = 2
singularity.enabled = true
singularity.cacheDir = "$HOME/.singularity/cache"

params {
    bam = "../data/alignments.bam"
    threads = 1
}
```
<br>

This tells nextflow how to run, and sets up inputs parameters we can use to supply values to the `SAMTOOLS_FLAGSTAT` process:

- The `bam` parameter is the input bam file we wish to analyse. 
- The `threads` parameter is an integer, and controls how many additional compute threads to use during runtime.

From here, we can refer to these inputs as `params.bam` / `params.threads` in other files.

<br>

**Creating Workflow & Passing Data** 

Now that we have the `nextflow.config` file set up, we will add a few lines to `samtools_flagstat.nf` to turn it into a workflow. 

Copy and paste the following lines at the top of `samtools_flagstat.nf`:

```
nextflow.enable.dsl=2

ch_bam = Channel.fromPath( params.bam )

workflow {
    SAMTOOLS_FLAGSTAT(
        ch_bam,             // input1 
        params.threads      // addthreads
    )
    SAMTOOLS_FLAGSTAT.out.output1.view()
}
```

The first line creates a nextflow `Channel` for our `bam` input. <br>
The `params.bam` global variable we set up previously is used to supply the 
path to our sample data.

The new `workflow {}` section declares the main workflow entry point. <br>
When we run this file, nextflow will look for this section and run the workflow contained within. 

In our case, the workflow only contains a single task, which runs the `SAMTOOLS_FLAGSTAT` process defined below the workflow section. We then supply input values to `SAMTOOLS_FLAGSTAT` using our `ch_bams` channel we created for `input1`, and `params.threads` for the `addthreads` input.

<br>

**Adding publishDir directive**

So that we can collect the output of `SAMTOOLS_FLAGSTAT` when it runs, we will add a `publishDir` directive to the process:

```
process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.13--h8c37831_0"
    publishDir "./outputs"
    
    ...

}
```

Now that we have set up `SAMTOOLS_FLAGSTAT` as a workflow, we can run it and check the output. 

<br>

**Running Our Workflow**

Ensure you are in the `translated/` working directory, where `nextflow.config` and `samtools_flagstat.nf` reside. 

```
cd translated
```

To run the workflow using our sample data, we can now write the following command: 
```
nextflow run samtools_flagstat.nf
```

Once completed, the check the `./outputs` folder inside `translated/`. 

If everything went well, you should see a single file called `output1.txt` with the following contents: 
```
200 + 0 in total (QC-passed reads + QC-failed reads)
200 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
25 + 0 mapped (12.50% : N/A)
25 + 0 primary mapped (12.50% : N/A)
200 + 0 paired in sequencing
100 + 0 read1
100 + 0 read2
0 + 0 properly paired (0.00% : N/A)
0 + 0 with itself and mate mapped
25 + 0 singletons (12.50% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

<br>

## Conclusion

In this tutorial we explored how to translate a simple Galaxy Tool to a Nextflow process. 

If needed, you can check the `final/` folder as a reference for the translated nextflow process and config. 

<br>






