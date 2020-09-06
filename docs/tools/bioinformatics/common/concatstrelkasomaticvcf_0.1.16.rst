:orphan:

Concat Strelka Somatic Vcf
====================================================

``ConcatStrelkaSomaticVcf`` · *0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.concat_strelkasomaticvcf import ConcatStrelkaSomaticVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "concatstrelkasomaticvcf_step",
           ConcatStrelkaSomaticVcf(
               headerVcfs=None,
               contentVcfs=None,
           )
       )
       wf.output("out", source=concatstrelkasomaticvcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for ConcatStrelkaSomaticVcf:

.. code-block:: bash

   # user inputs
   janis inputs ConcatStrelkaSomaticVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       contentVcfs:
       - contentVcfs_0.vcf.gz
       - contentVcfs_1.vcf.gz
       headerVcfs:
       - headerVcfs_0.vcf.gz
       - headerVcfs_1.vcf.gz




5. Run ConcatStrelkaSomaticVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       ConcatStrelkaSomaticVcf





Information
------------

:ID: ``ConcatStrelkaSomaticVcf``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.16
:Container: biocontainers/vcftools:v0.1.16-1-deb_cv1
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ===========================  ========  ==========  ===============
name            type                         prefix      position  documentation
==============  ===========================  ========  ==========  ===============
headerVcfs      Array<CompressedIndexedVCF>                     1
contentVcfs     Array<CompressedIndexedVCF>                     4
outputFilename  Optional<Filename>           >                  6
==============  ===========================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task ConcatStrelkaSomaticVcf {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] headerVcfs
       Array[File] headerVcfs_tbi
       Array[File] contentVcfs
       Array[File] contentVcfs_tbi
       String? outputFilename
     }
     command <<<
       set -e
        \
         vcf-merge \
         ~{"'" + sep("' '", headerVcfs) + "'"} \
         | grep '^##' > header.vcf; \
         vcf-concat \
         ~{"'" + sep("' '", contentVcfs) + "'"} \
         | grep -v '^##' > content.vcf; cat header.vcf content.vcf \
         > ~{select_first([outputFilename, "generated.strelka.vcf"])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/vcftools:v0.1.16-1-deb_cv1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.strelka.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.0
   label: Concat Strelka Somatic Vcf

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/vcftools:v0.1.16-1-deb_cv1

   inputs:
   - id: headerVcfs
     label: headerVcfs
     type:
       type: array
       items: File
     inputBinding:
       position: 1
   - id: contentVcfs
     label: contentVcfs
     type:
       type: array
       items: File
     inputBinding:
       position: 4
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.strelka.vcf
     inputBinding:
       prefix: '>'
       position: 6
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.strelka.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: vcf-merge
     shellQuote: false
   - position: 2
     valueFrom: "| grep '^##' > header.vcf;"
     shellQuote: false
   - position: 3
     valueFrom: vcf-concat
     shellQuote: false
   - position: 5
     valueFrom: "| grep -v '^##' > content.vcf; cat header.vcf content.vcf"
     shellQuote: false
   id: ConcatStrelkaSomaticVcf


