:orphan:

Filter Vardict Somatic Vcf
====================================================

``FilterVardictSomaticVcf`` · *2 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.filtervardictsomaticvcf import FilterVardictSomaticVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "filtervardictsomaticvcf_step",
           FilterVardictSomaticVcf(

           )
       )
       wf.output("out", source=filtervardictsomaticvcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FilterVardictSomaticVcf:

.. code-block:: bash

   # user inputs
   janis inputs FilterVardictSomaticVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run FilterVardictSomaticVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FilterVardictSomaticVcf





Information
------------

:ID: ``FilterVardictSomaticVcf``
:URL: *No URL to the documentation was provided*
:Versions: v1.9
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
:Authors: Jiaan Yu, Michael Franklin
:Citations: None
:Created: 2020-06-04
:Updated: 2020-11-09


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
vcf             Optional<VCF>                          1
outputFilename  Optional<Filename>  -o                 3
==============  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task FilterVardictSomaticVcf {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? vcf
       String? outputFilename
     }
     command <<<
       set -e
        \
         bcftools filter -e 'STATUS="GERMLINE"' -o - \
         ~{if defined(vcf) then ("'" + vcf + "'") else ""} \
         | bcftools filter -i 'FILTER=="PASS"' \
         -o ~{select_first([outputFilename, "~{basename(vcf, ".vcf")}.filter.vcf"])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(vcf, ".vcf")}.filter.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Filter Vardict Somatic Vcf
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/bcftools:v1.9-1-deb_cv1

   inputs:
   - id: vcf
     label: vcf
     type:
     - File
     - 'null'
     inputBinding:
       position: 1
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.filter.vcf
     inputBinding:
       prefix: -o
       position: 3
       valueFrom: |-
         $(inputs.vcf ? inputs.vcf.basename.replace(/.vcf$/, "") : "generated").filter.vcf
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: |-
         $(inputs.vcf ? inputs.vcf.basename.replace(/.vcf$/, "") : "generated").filter.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: "bcftools filter -e 'STATUS=\"GERMLINE\"' -o - "
     shellQuote: false
   - position: 2
     valueFrom: "| bcftools filter -i 'FILTER==\"PASS\"'"
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: FilterVardictSomaticVcf


