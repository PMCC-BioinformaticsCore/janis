:orphan:

BCFTools: Sort
=============================

``bcftoolssort`` · *1 contributor · 1 version*

About:   Sort VCF/BCF file.
Usage:   bcftools sort [OPTIONS] <FILE.vcf>


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.sort.versions import BcfToolsSort_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolssort_step",
           BcfToolsSort_1_9(
               vcf=None,
           )
       )
       wf.output("out", source=bcftoolssort_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolssort:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolssort > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: null




5. Run bcftoolssort with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolssort





Information
------------

:ID: ``bcftoolssort``
:URL: *No URL to the documentation was provided*
:Versions: v1.9
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-09
:Updated: 2019-07-11


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     Gzipped<VCF>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

==============  ========================  =============  ==========  =======================================================================================
name            type                      prefix           position  documentation
==============  ========================  =============  ==========  =======================================================================================
vcf             Union<VCF, Gzipped<VCF>>                          1  The VCF file to sort
outputFilename  Optional<Filename>        --output-file              (-o) output file name [stdout]
outputType      Optional<String>          --output-type              (-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
tempDir         Optional<String>          --temp-dir                 (-T) temporary files [/tmp/bcftools-sort.XXXXXX/]
==============  ========================  =============  ==========  =======================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcftoolssort {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? outputFilename
       String? outputType
       String? tempDir
     }
     command <<<
       set -e
       bcftools sort \
         --output-file '~{select_first([outputFilename, "generated.sorted.vcf.gz"])}' \
         ~{if defined(select_first([outputType, "z"])) then ("--output-type '" + select_first([outputType, "z"]) + "'") else ""} \
         ~{if defined(tempDir) then ("--temp-dir '" + tempDir + "'") else ""} \
         ~{vcf}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.sorted.vcf.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BCFTools: Sort'
   doc: "About:   Sort VCF/BCF file.\nUsage:   bcftools sort [OPTIONS] <FILE.vcf>"

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/bcftools:v1.9-1-deb_cv1

   inputs:
   - id: vcf
     label: vcf
     doc: The VCF file to sort
     type: File
     inputBinding:
       position: 1
   - id: outputFilename
     label: outputFilename
     doc: (-o) output file name [stdout]
     type:
     - string
     - 'null'
     default: generated.sorted.vcf.gz
     inputBinding:
       prefix: --output-file
   - id: outputType
     label: outputType
     doc: |-
       (-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
     type: string
     default: z
     inputBinding:
       prefix: --output-type
   - id: tempDir
     label: tempDir
     doc: (-T) temporary files [/tmp/bcftools-sort.XXXXXX/]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --temp-dir

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.sorted.vcf.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bcftools
   - sort
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcftoolssort


