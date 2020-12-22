:orphan:

CNVKit
======

``CNVKit`` · *1 contributor · 1 version*


        A command-line toolkit and Python library for detecting copy number variants 
        and alterations genome-wide from high-throughput sequencing.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.ucsf.cnvkit.cnvkit_0_9_6 import CNVKit_0_9_6

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cnvkit_step",
           CNVKit_0_9_6(
               reference=None,
           )
       )

    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for CNVKit:

.. code-block:: bash

   # user inputs
   janis inputs CNVKit > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference




5. Run CNVKit with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       CNVKit





Information
------------

:ID: ``CNVKit``
:URL: `https://github.com/etal/cnvkit <https://github.com/etal/cnvkit>`_
:Versions: 0.9.6
:Container: etal/cnvkit:0.9.6
:Authors: Michael Franklin
:Citations: Talevich, E., Shain, A.H., Botton, T., & Bastian, B.C. (2014). CNVkit: Genome-wide copy number detection and visualization from targeted sequencing. PLOS Computational Biology 12(4):e1004873
:DOI: 10.1371/journal.pcbi.1004873
:Created: 2019-07-03 00:00:00
:Updated: 2019-07-03 00:00:00


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
======  ======  ===============


Additional configuration (inputs)
---------------------------------

===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
name             type                prefix               position    documentation
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================
reference        File                --reference                      REFERENCE Copy number reference file (.cnn).
outputDirectory  Optional<Filename>  --output-dir                     DIRECTORY Output directory.
method           Optional<String>    --method                         (-m) {hybrid,amplicon,wgs} Sequencing protocol: hybridization capture ('hybrid'), targeted amplicon sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines whether and how to use antitarget bins. [Default: hybrid]
maleReference    Optional<String>    --male-reference                 (-y, --haploid-x-reference) Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).
countReads       Optional<String>    --count-reads                    (-c) Get read depths by counting read midpoints within each bin. (An alternative algorithm).
dropLowCoverage  Optional<String>    --drop-low-coverage              Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
processes        Optional<String>    --processes                      (-p) [PROCESSES] Number of subprocesses used to running each of the BAM files in parallel. Without an argument, use the maximum number of available CPUs. [Default: process each BAM in serial]
rscriptPath      Optional<String>    --rscript-path                   Path to the Rscript excecutable to use for running R code. Use this option to specify a non-default R installation. [Default: Rscript]
===============  ==================  ===================  ==========  =====================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task CNVKit {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? outputDirectory
       File reference
       String? method
       String? maleReference
       String? countReads
       String? dropLowCoverage
       String? processes
       String? rscriptPath
     }
     command <<<
       set -e
       cnvkit.py batch \
         --output-dir '~{select_first([outputDirectory, "generated"])}' \
         --reference '~{reference}' \
         ~{if defined(method) then ("--method '" + method + "'") else ""} \
         ~{if defined(maleReference) then ("--male-reference '" + maleReference + "'") else ""} \
         ~{if defined(countReads) then ("--count-reads '" + countReads + "'") else ""} \
         ~{if defined(dropLowCoverage) then ("--drop-low-coverage '" + dropLowCoverage + "'") else ""} \
         ~{if defined(processes) then ("--processes '" + processes + "'") else ""} \
         ~{if defined(rscriptPath) then ("--rscript-path '" + rscriptPath + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "etal/cnvkit:0.9.6"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: CNVKit
   doc: |2-

             A command-line toolkit and Python library for detecting copy number variants 
             and alterations genome-wide from high-throughput sequencing.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: etal/cnvkit:0.9.6

   inputs:
   - id: outputDirectory
     label: outputDirectory
     doc: DIRECTORY Output directory.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --output-dir
   - id: reference
     label: reference
     doc: REFERENCE Copy number reference file (.cnn).
     type: File
     inputBinding:
       prefix: --reference
   - id: method
     label: method
     doc: |-
       (-m) {hybrid,amplicon,wgs} Sequencing protocol: hybridization capture ('hybrid'), targeted amplicon sequencing ('amplicon'), or whole genome sequencing ('wgs'). Determines whether and how to use antitarget bins. [Default: hybrid]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --method
   - id: maleReference
     label: maleReference
     doc: |-
       (-y, --haploid-x-reference) Use or assume a male reference (i.e. female samples will have +1 log-CNR of chrX; otherwise male samples would have -1 chrX).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --male-reference
   - id: countReads
     label: countReads
     doc: |2-
        (-c) Get read depths by counting read midpoints within each bin. (An alternative algorithm).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --count-reads
   - id: dropLowCoverage
     label: dropLowCoverage
     doc: |-
       Drop very-low-coverage bins before segmentation to avoid false-positive deletions in poor-quality tumor samples.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --drop-low-coverage
   - id: processes
     label: processes
     doc: |-
       (-p) [PROCESSES] Number of subprocesses used to running each of the BAM files in parallel. Without an argument, use the maximum number of available CPUs. [Default: process each BAM in serial]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --processes
   - id: rscriptPath
     label: rscriptPath
     doc: |-
       Path to the Rscript excecutable to use for running R code. Use this option to specify a non-default R installation. [Default: Rscript]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --rscript-path

   outputs: []
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - cnvkit.py
   - batch
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: CNVKit


