# Tutorial 3: CWL Tool -> Nextflow

<br>

## Sections

- [Introduction](#introduction)
- [Running Janis Translate](#running-janis-translate)
- [Manual Adjustments to Translated Tool](#manual-adjustments-to-translated-tool)
- [Running Translated Tool as a Workflow](#running-translated-tool-as-a-workflow)
- [Conclusion](#conclusion)

<br>

## Introduction

This tutorial demonstrates translation of a `gatk HaplotypeCaller` tool from CWL to Nextflow using `janis translate`. <br>

The source CWL file used in this section is taken from the [McDonnell Genome Institute](https://www.genome.wustl.edu/) (MGI) [analysis-workflows](https://github.com/genome/analysis-workflows) repository. 

This repository stores publically available analysis pipelines for genomics data. <br>
It is a fantastic piece of research software, and the authors thank MGI for their contribution to open-source research software. 

The software tool encapsulated by the source CWL - [gatk_haplotype_caller.cwl](https://github.com/genome/analysis-workflows/blob/master/definitions/tools/gatk_haplotype_caller.cwl) - displays summary information for an alignment file. 

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
wget https://zenodo.org/record/8275567/files/tutorial3.tar.gz
tar -xvf tutorial3.tar.gz
```

After you have decompressed the tar archive, change into the new directory: 

```
cd tutorial3
```

Inside this folder we have the following structure: 
```
tutorial3
├── data
│   ├── 2895499223_sorted.bam
│   ├── 2895499223_sorted.bam.bai
│   ├── chr17_test.dict
│   ├── chr17_test.fa
│   ├── chr17_test.fa.fai
│   ├── chr17_test_dbsnp.vcf.gz
│   └── chr17_test_dbsnp.vcf.gz.tbi
├── final
│   ├── gatk_haplotype_caller.nf
│   └── nextflow.config
└── source
    └── gatk_haplotype_caller.cwl
``` 

We will translate the `source/gatk_haplotype_caller.cwl` tool into nextflow using janis, then will test out our translation using the files inside the indexed bam file in the `data/` folder. 


<br>

## Running Janis Translate

To translate a tool / workflow,  we use `janis translate`.

```
janis translate --from <src> --to <dest> <filepath>
```

The `--from` argument specifies the workflow language of the source file(s), and `--to` specifies the destination we want to translate to. 

In our case, this will be `--from cwl --to nextflow`.

The `<filepath>` argument is the source file we will translate. 

In this tutorial, the filepath will be `source/gatk_haplotype_caller.cwl`

<br>

To translate our cwl tool, run the following command:
```
singularity exec ~/janis.sif janis translate --from cwl --to nextflow source/gatk_haplotype_caller.cwl
```

You will see a folder called `translated/` appear, and a nextflow process called `gatk_haplotype_caller.nf` will be present inside. 

<br>

## Manual Adjustments to Translated Tool

The `translated/gatk_haplotype_caller.nf` file should be similar to the following: 

```
nextflow.enable.dsl = 2

process GATK_HAPLOTYPE_CALLER {
    
    container "broadinstitute/gatk:4.1.8.1"

    input:
    path bam
    path reference
    path dbsnp_vcf, stageAs: 'dbsnp_vcf/*'
    val intervals
    val gvcf_gq_bands
    val emit_reference_confidence
    val contamination_fraction
    val max_alternate_alleles
    val ploidy
    val read_filter

    output:
    tuple path("output.g.vcf.gz"), path("output.g.vcf.gz.tbi"), emit: gvcf

    script:
    def bam = bam[0]
    def dbsnp_vcf = dbsnp_vcf[0] != null ? "--dbsnp ${dbsnp_vcf[0]}" : ""
    def reference = reference[0]
    def gvcf_gq_bands_joined = gvcf_gq_bands != params.NULL_VALUE ? "-GQB " + gvcf_gq_bands.join(' ') : ""
    def intervals_joined = intervals.join(' ')
    def contamination_fraction = contamination_fraction != params.NULL_VALUE ? "-contamination ${contamination_fraction}" : ""
    def max_alternate_alleles = max_alternate_alleles != params.NULL_VALUE ? "--max_alternate_alleles ${max_alternate_alleles}" : ""
    def ploidy = ploidy != params.NULL_VALUE ? "-ploidy ${ploidy}" : ""
    def read_filter = read_filter != params.NULL_VALUE ? "--read_filter ${read_filter}" : ""
    """
    /gatk/gatk --java-options -Xmx16g HaplotypeCaller \
    -R ${reference} \
    -I ${bam} \
    -ERC ${emit_reference_confidence} \
    ${gvcf_gq_bands_joined} \
    -L ${intervals_joined} \
    ${dbsnp_vcf} \
    ${contamination_fraction} \
    ${max_alternate_alleles} \
    ${ploidy} \
    ${read_filter} \
    -O "output.g.vcf.gz"
    """

}
```

We can see that this nextflow process has a multiple inputs, single output, and calls `gatk HaplotypeCaller` using the input data we supply to the process.  

<br>

> Notes on translation: <br><br>
> **(1)** `def bam = bam[0]`<br><br>
> This pattern is used in the script block to handle datatypes with secondary files. <br>
> The `bam` input is an indexed bam type, so requires a `.bai` file to also be present in the working directory. <br><br>
> To facilitate this, we supply the `bam` input as a tuple of 2 files:<br> 
> `[filename.bam, filename.bam.bai]`. <br><br>
> `def bam = bam[0]` is used so when we reference `${bam}` in the script body, we are refering to the `.bam` file in that list. <br><br>
> **(2)** `def ploidy = ploidy != params.NULL_VALUE ? "-ploidy ${ploidy}" : ""` <br><br>
> This is how we handle optional val inputs in nextflow. <br>
> Nextflow doesn't like ***null*** values to be passed to process inputs, so we use 2 different tricks to make ***optional*** inputs possible.<br><br>
> For `val` inputs we set up a `NULL_VALUE` param in `nextflow.config` which we use as a placeholder.  <br>
> For `path` inputs (ie files and directories) we set up a ***null*** file rather than passing ***null*** directly. <br>
> This ensures that the file is staged correctly in the working directory when an actual filepath is provided. <br><br>
> The format above is a ternary operator of the form `def myvar = cond_check ? cond_true : cond_false`.<br>
> We are redefining the `ploidy` string variable so that it will be correctly formatted when used in the script block.<br>
> If we supply an actual value, it will take the form `"-ploidy ${ploidy}"`. <br>
> If we supply the placeholder `params.NULL_VALUE` value, it will be a blank string `""`. <br><br>
> **(3)** `def intervals_joined = intervals.join(' ')` <br><br>
> Templates our `intervals` list of strings to a single string.<br>
> Each item is joined by a space: eg `["hello", "there"]` -> `"hello there"`. <br><br>
> **(4)** `def dbsnp_vcf = dbsnp_vcf[0] != null ? "--dbsnp ${dbsnp_vcf[0]}" : ""`<br><br>
> Same as **(2)** except uses a different check for null value, and selects the primary `.vcf.gz` file. <br><br>
> As the `path dbsnp_vcf` input is optional and consists of a `.vcf.gz` and `.vcf.gz.tbi` file, `dbsnp_vcf[0] != null` checks if the primary `.vcf.gz` file is present. <br><br>
> If true, it templates a string we can use in the script body to provide the required argument: `--dbsnp ${dbsnp_vcf[0]}`.<br>
> If false, it templates an empty string.<br><br>
> This ensures that when we use `${dbsnp_vcf}` in the script body, it will be formatted correctly for either case.<br><br>

<br>

This translation is correct for the `gatk_haplotype_caller.cwl` file and needs no adjusting. <br>
Have a look at the source CWL file to see how they match up. 

<br>

## Running Translated Tool as a Workflow

**Collecting Process Outputs**

Let's add a `publishDir` directive to our translated process so that we can capture the outputs of this process.

```
process GATK_HAPLOTYPE_CALLER {

    container "broadinstitute/gatk:4.1.8.1"
    publishDir "./outputs"    
    ...
}                               
```

Nextflow allows us to capture the outputs created by a process using the `publishDir` directive seen above. 

<br>

>NOTE<br>
> Our `gatk_haplotype_caller.nf` has a container directive with the form `container "broadinstitute/gatk:4.1.8.1"` provided, so Nextflow will handle this singularity image download and will use the specified image when running this process (provided that `singularity.enabled` is set to `true` in the nextflow config). <br>
> If you'd like to use an already available container, you can modify this container directive to `container "/cvmfs/singularity.galaxyproject.org/all/gatk4:4.1.8.1--py38_0"`.

<br>

**Setting up nextflow.config**

To run this process, we will set up a `nextflow.config` file and add some lines to the top of our process definition to turn it into a workflow.

Create a new file called `nextflow.config` in the `translated/` folder alongside `gatk_haplotype_caller.nf`. 

Copy and paste the following code into your `nextflow.config` file: 

```
nextflow.enable.dsl = 2
singularity.enabled = true
singularity.cacheDir = "$HOME/.singularity/cache"

params {
    NULL_VALUE = 'NULL'

    bam = [
        '../data/2895499223_sorted.bam',
        '../data/2895499223_sorted.bam.bai',
    ]
    dbsnp = [
        '../data/chr17_test_dbsnp.vcf.gz',
        '../data/chr17_test_dbsnp.vcf.gz.tbi',
    ]
    reference = [
        '../data/chr17_test.fa',
        '../data/chr17_test.dict',
        '../data/chr17_test.fa.fai',
    ]
    gvcf_gq_bands = NULL_VALUE
    intervals = ["chr17"]
    emit_reference_confidence = 'GVCF'
    contamination_fraction = NULL_VALUE
    max_alternate_alleles = NULL_VALUE
    ploidy = NULL_VALUE
    read_filter = NULL_VALUE
}
```
This tells nextflow how to run, and sets up the sample data as inputs.

*file inputs*

The `bam` parameter is a list which provides paths to the `.bam` and `.bai` sample data we will use to test the nextflow translation. From here, we can refer to the indexed bam input as `params.bam` in other files. The `dbsnp` and `reference` params follow this same pattern. 

*non-file inputs*

We also set up a `NULL_VALUE` param which we use as a *placeholder* for a null value. <br> 
In this case we are providing null values for the `gvcf_gq_bands`, `contamination_fraction`, `max_alternate_alleles`, `ploidy` and `read_filter` inputs as they are all optional.

<br>

**Creating Workflow & Passing Data** 

Now that we have the `nextflow.config` file set up, we will add a few lines to `gatk_haplotype_caller.nf` to turn it into a workflow. 

Copy and paste the following lines at the top of `gatk_haplotype_caller.nf`:

```
ch_bam = Channel.fromPath( params.bam ).toList()
ch_dbsnp = Channel.fromPath( params.dbsnp ).toList()
ch_reference = Channel.fromPath( params.reference ).toList()

workflow {

    GATK_HAPLOTYPE_CALLER(
        ch_bam,
        ch_reference,
        ch_dbsnp,
        params.intervals,
        params.gvcf_gq_bands,
        params.emit_reference_confidence,
        params.contamination_fraction,
        params.max_alternate_alleles,
        params.ploidy,
        params.read_filter,
    )

}
```

The first 3 lines create nextflow `Channels` for our `bam`, `dbsnp`, and `reference` inputs and ensures they are lists.

The `Channel.toList()` aspect collects our files into a list, as the primary & secondary files for these datatypes must be passed together as a tuple.

The `params.bam`, `params.dbsnp` and `params.reference` global variables we set up previously are used to supply the paths to our sample data for these channels.

The new `workflow {}` section declares the main workflow entry point. <br>
When we run this file, nextflow will look for this section and run the workflow contained within. 

In our case, the workflow only contains a single task, which runs the `GATK_HAPLOTYPE_CALLER` process defined below the workflow section. We call `GATK_HAPLOTYPE_CALLER` by feeding inputs in the correct order, using the channels we declared at the top of the file, and variables we set up in the global `params` object. 

<br>

**Running Our Workflow**

Ensure you are in the `translated/` working directory, where `nextflow.config` and `gatk_haplotype_caller.nf` reside. 
```
cd translated
```

To run the workflow using our sample data, we can now write the following command: 
```
nextflow run gatk_haplotype_caller.nf
```

Nextflow will automatically check if there is a `nextflow.config` file in the working directory, and if so will use that to configure itself. Our inputs are supplied in `nextflow.config` alongside the dsl2 & singularity config, so it should run without issue. 

Once completed, we can check the `outputs/` folder to view our results. <br>
If everything went well, the `outputs/` folder should contain 2 files: 

- `output.g.vcf.gz`
- `output.g.vcf.gz.tbi`

<br>

## Conclusion

In this section we explored how to translate the `gatk_haplotype_caller` CWL CommandLineTool to a Nextflow process. 

This is a more real-world situation, where the CommandLineTool has multiple inputs, secondary files, and optionality. 

If needed, you can check the `tutorial3/final` folder which contains the nextflow files we created in this tutorial as reference.  

<br>