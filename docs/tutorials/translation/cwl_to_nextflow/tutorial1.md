# Tutorial 1: Simple CWL Tool -> Nextflow

<br>

## Sections

- [Introduction](#introduction)
- [Running Janis Translate](#running-janis-translate)
- [Manual Adjustments](#manual-adjustments-to-translated-tool)
- [Running Translated Tool](#running-translated-tool-as-a-workflow)
- [Conclusion](#conclusion)

<br>

## Introduction

This section demonstrates translation of a basic `samtools flagstat` tool from CWL to Nextflow using `janis translate`. 

The CWL tool used in this section - [samtools_flagstat.cwl](https://github.com/genome/analysis-workflows/blob/master/definitions/tools/samtools_flagstat.cwl) -  is taken from the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) [analysis-workflows](https://github.com/genome/analysis-workflows) repository. 

This resource stores publically available analysis pipelines for genomics data. <br>
It is a fantastic piece of research software, and the authors thank MGI for their contribution to open-source research software. 

The underlying software run by this tool - [Samtools Flagstat](http://www.htslib.org/doc/samtools-flagstat.html) - displays summary information for an alignment file.

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

This archive contains the source CWL file, sample data, and finished translations as a reference.

Run the following commands to download & decompress the zenodo archive:
```
wget https://zenodo.org/record/8275567/files/tutorial1.tar.gz
tar -xvf tutorial1.tar.gz
```

After you have decompressed the tar archive, change into the new directory: 

```
cd tutorial1
```

Inside this folder we have the following structure: 
```
tutorial1
├── data
│   ├── 2895499223_sorted.bam
│   └── 2895499223_sorted.bam.bai
├── final
│   ├── nextflow.config
│   └── samtools_flagstat.nf
└── source
    └── samtools_flagstat.cwl
``` 

We will translate the `source/samtools_flagstat.cwl` tool into nextflow using janis, then will test out our translation using the indexed bam file in the `data/` folder. 

<br>

## Running Janis Translate

To translate a tool / workflow,  we use `janis translate`.

```
janis translate --from <src> --to <dest> <filepath>
```

The `--from` argument specifies the workflow language of the source file(s), and `--to` specifies the destination we want to translate to. 

In our case, this will be `--from cwl --to nextflow`.

The `<filepath>` argument is the source file we will translate. 

In this tutorial, the filepath will be `source/samtools_flagstat.cwl`

<br>

To translate our cwl tool, run the following command:
```
singularity exec ~/janis.sif janis translate --from cwl --to nextflow source/samtools_flagstat.cwl
```

You will see a folder called `translated/` appear, and a nextflow process called `samtools_flagstat.nf` will be present inside. 

<br>

## Manual Adjustments to Translated Tool

The `translated/samtools_flagstat.nf` file should be similar to the following: 

```
nextflow.enable.dsl = 2

process SAMTOOLS_FLAGSTAT {
    
    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"

    input:
    path bam

    output:
    path "${bam[0]}.flagstat", emit: flagstats

    script:
    def bam = bam[0]
    """
    /usr/local/bin/samtools flagstat \
    ${bam} \
    > ${bam}.flagstat \
    """

}
```

We can see that this nextflow process has a single input, a single output, and calls `samtools flagstat` on the input `bam`. 

We can also see that a container image is available for this tool. In the next section we will run this process using some sample data and the specified container. 

This translation is correct for the `samtools_flagstat.cwl` file and needs no adjusting. <br>
Have a look at the source CWL file to see how they match up. 

<br>

> Note: <br>
> `def bam = bam[0]` in the script block is used to handle datatypes with secondary files. <br>
> The `bam` input is an indexed bam type, so requires a `.bai` file to also be present in the working directory alongside the `.bam` file. <br><br>
> For this reason the `bam` input is supplied as an Array with 2 files - the `.bam` and the `.bai`. <br>
> Here the `def bam = bam[0]` is used so that `bam` refers to the `.bam` file in that Array. 

<br>

## Running Translated Tool as a Workflow

**Collecting Process Outputs**

Let's add a `publishDir` directive to our translated process so that we can capture the outputs of this process.

```
process SAMTOOLS_FLAGSTAT {

    container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"
    publishDir "./outputs"                                          <-
    ....

}                          
```

Nextflow allows us to capture the outputs created by a process using the `publishDir` directive seen above. 

<br>

**Setting up nextflow.config**

To run this process, we will set up a `nextflow.config` file and add some lines to the top of our process definition to turn it into a workflow.

Create a new file called `nextflow.config` in the `translated/` folder alongside `samtools_flagstat.nf`. 

Copy and paste the following code into your `nextflow.config` file: 

```
nextflow.enable.dsl = 2
singularity.enabled = true
singularity.cacheDir = "$HOME/.singularity/cache"

params {

    bam = [
        '../data/2895499223_sorted.bam',
        '../data/2895499223_sorted.bam.bai',
    ]

}
```

This tells nextflow how to run, and sets up an input parameter for our indexed bam input.

The `bam` parameter is a list which provides paths to the `.bam` and `.bai` sample data we will use to test the nextflow translation. From here, we can refer to the indexed bam input as `params.bam` in other files.

<br>

>NOTE<br>
>`nextflow.enable.dsl = 2` ensures that we are using the dsl2 nextflow syntax which is the current standard. <br><br>
>`singularity.enabled = true` tells nextflow to run processes using singularity. Our `samtools_flagstat.nf` has a directive with the form `container "quay.io/biocontainers/samtools:1.11--h6270b1f_0"` provided, so it will use the specified image when running this process. <br> <br>
>`singularity.cacheDir = "$HOME/.singularity/cache"` tells nextflow where singularity images are stored. <br>
> Nextflow will handle the singularity image download and stored it in the cache specified above.

<br>

**Creating Workflow & Passing Data** 

Now that we have the `nextflow.config` file set up, we will add a few lines to `samtools_flagstat.nf` to turn it into a workflow. 

Copy and paste the following lines at the top of `samtools_flagstat.nf`:

```
ch_bam = Channel.fromPath( params.bam ).toList()

workflow {
    SAMTOOLS_FLAGSTAT(ch_bam)
}
```

The first line creates a nextflow `Channel` for our `bam` input and ensures it is a list. <br>
The `Channel.toList()` part collects our files into a list, as both the `.bam` and `.bai` files must be passed together. <br>
The `params.bam` global variable we set up previously is used to supply the paths to our sample data.

The new `workflow {}` section declares the main workflow entry point. <br>
When we run this file, nextflow will look for this section and run the workflow contained within. 

In our case, the workflow only contains a single task, which runs the `SAMTOOLS_FLAGSTAT` process defined below the workflow section. The single `SAMTOOLS_FLAGSTAT` input is being passed data from our `ch_bam` channel we declared at the top of the file. 

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

Nextflow will automatically check if there is a `nextflow.config` file in the working directory, and if so will use that to configure itself. Our inputs are supplied in `nextflow.config` alongside the dsl2 & singularity config, so it should run without issue. 

Once completed, we can check the `outputs/` folder to view our results. 

If everything went well, there should be a file called `2895499223_sorted.bam.flagstat` with the following contents:

```
2495 + 0 in total (QC-passed reads + QC-failed reads)
0 + 0 secondary
1 + 0 supplementary
0 + 0 duplicates
2480 + 0 mapped (99.40% : N/A)
2494 + 0 paired in sequencing
1247 + 0 read1
1247 + 0 read2
2460 + 0 properly paired (98.64% : N/A)
2464 + 0 with itself and mate mapped
15 + 0 singletons (0.60% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)

```

<br>

## Conclusion

In this tutorial we explored how to translate a simple CWL tool to a Nextflow process. 

If you ran into difficulty, you can check the `tutorial1/final/` folder. <br>
This folder contains the nextflow files we created in this tutorial as a reference. 

<br>