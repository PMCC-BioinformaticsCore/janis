:orphan:

Gene Coverage Per Sample
================================================

``geneCoveragePerSample`` · *1 contributor · 2 versions*

usage: gene_coverage_per_sample.py [-h] [-l LIST] [-n NAME] [-p PATH] [-b BED]
                                   [-g GENE] [-r REGION] [-f FOLDS] [-d]
                                   [-t THREADS]

Gene or region coverage of bam

optional arguments:
  -h, --help            show this help message and exit
  -l LIST, --list LIST  List file: A tsv file contains SampleName
                        PathToBedtoolsOutput on each line
  -n NAME, --name NAME  Sample name if list not used
  -p PATH, --path PATH  Path to bedtools output if list not used
  -b BED, --bed BED     (Deprecated option) Bed file
  -g GENE, --gene GENE  Output gene file
  -r REGION, --region REGION
                        Output region file
  -f FOLDS, --folds FOLDS
                        Folds, quoted and commna sepparated, default
                        1,10,20,100
  -d, --remove_duplicates
                        (Deprecated option) Remove marked duplicates in
                        analysis, default:false
  -t THREADS, --threads THREADS
                        number of threads, default:32
        


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.genecovpersample.versions import GeneCoveragePerSample_0_0_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "genecoveragepersample_step",
           GeneCoveragePerSample_0_0_7(

           )
       )
       wf.output("geneFileOut", source=genecoveragepersample_step.geneFileOut)
       wf.output("regionFileOut", source=genecoveragepersample_step.regionFileOut)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for geneCoveragePerSample:

.. code-block:: bash

   # user inputs
   janis inputs geneCoveragePerSample > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run geneCoveragePerSample with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       geneCoveragePerSample





Information
------------

:ID: ``geneCoveragePerSample``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/performance>`_
:Versions: 0.0.8, 0.0.7
:Container: michaelfranklin/pmacutil:0.0.7
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-04-03 00:00:00
:Updated: 2020-04-03 00:00:00


Outputs
-----------

=============  ========  ===============
name           type      documentation
=============  ========  ===============
geneFileOut    TextFile
regionFileOut  TextFile
=============  ========  ===============


Additional configuration (inputs)
---------------------------------

==================  ==================  =========  ==========  ========================================================
name                type                prefix     position    documentation
==================  ==================  =========  ==========  ========================================================
listFile            Optional<File>      --list                 List file: A tsv file contains SampleName	PathToBedtoolsOutput on each line
sampleName          Optional<String>    --name                 Sample name if list not used
bedtoolsOutputPath  Optional<File>      --path                 Path to bedtools output if list not used
outputGeneFile      Optional<Filename>  --gene                 Output gene file
outputRegionFile    Optional<Filename>  --region               Output region file
fold                Optional<String>    --fold                 Folds, quoted and commna sepparated, default 1,10,20,100
threads             Optional<Integer>   --threads              number of threads, default:32
==================  ==================  =========  ==========  ========================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task geneCoveragePerSample {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? listFile
       String? sampleName
       File? bedtoolsOutputPath
       String? outputGeneFile
       String? outputRegionFile
       String? fold
       Int? threads
     }
     command <<<
       set -e
       gene_coverage_per_sample.py \
         ~{if defined(listFile) then ("--list '" + listFile + "'") else ""} \
         ~{if defined(sampleName) then ("--name '" + sampleName + "'") else ""} \
         ~{if defined(bedtoolsOutputPath) then ("--path '" + bedtoolsOutputPath + "'") else ""} \
         --gene '~{select_first([outputGeneFile, "generated.gene.txt"])}' \
         --region '~{select_first([outputRegionFile, "generated.region.txt"])}' \
         ~{if defined(fold) then ("--fold '" + fold + "'") else ""} \
         ~{if defined(threads) then ("--threads " + threads) else ''}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.0.7"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File geneFileOut = select_first([outputGeneFile, "generated.gene.txt"])
       File regionFileOut = select_first([outputRegionFile, "generated.region.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Gene Coverage Per Sample
   doc: |-
     usage: gene_coverage_per_sample.py [-h] [-l LIST] [-n NAME] [-p PATH] [-b BED]
                                        [-g GENE] [-r REGION] [-f FOLDS] [-d]
                                        [-t THREADS]

     Gene or region coverage of bam

     optional arguments:
       -h, --help            show this help message and exit
       -l LIST, --list LIST  List file: A tsv file contains SampleName
                             PathToBedtoolsOutput on each line
       -n NAME, --name NAME  Sample name if list not used
       -p PATH, --path PATH  Path to bedtools output if list not used
       -b BED, --bed BED     (Deprecated option) Bed file
       -g GENE, --gene GENE  Output gene file
       -r REGION, --region REGION
                             Output region file
       -f FOLDS, --folds FOLDS
                             Folds, quoted and commna sepparated, default
                             1,10,20,100
       -d, --remove_duplicates
                             (Deprecated option) Remove marked duplicates in
                             analysis, default:false
       -t THREADS, --threads THREADS
                             number of threads, default:32
          

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.7

   inputs:
   - id: listFile
     label: listFile
     doc: "List file: A tsv file contains SampleName\tPathToBedtoolsOutput on each line"
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --list
   - id: sampleName
     label: sampleName
     doc: Sample name if list not used
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --name
   - id: bedtoolsOutputPath
     label: bedtoolsOutputPath
     doc: Path to bedtools output if list not used
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --path
   - id: outputGeneFile
     label: outputGeneFile
     doc: Output gene file
     type:
     - string
     - 'null'
     default: generated.gene.txt
     inputBinding:
       prefix: --gene
   - id: outputRegionFile
     label: outputRegionFile
     doc: Output region file
     type:
     - string
     - 'null'
     default: generated.region.txt
     inputBinding:
       prefix: --region
   - id: fold
     label: fold
     doc: Folds, quoted and commna sepparated, default 1,10,20,100
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --fold
   - id: threads
     label: threads
     doc: number of threads, default:32
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads

   outputs:
   - id: geneFileOut
     label: geneFileOut
     type: File
     outputBinding:
       glob: generated.gene.txt
       loadContents: false
   - id: regionFileOut
     label: regionFileOut
     type: File
     outputBinding:
       glob: generated.region.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: gene_coverage_per_sample.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: geneCoveragePerSample


