:orphan:

VcfLib: VcfFixUp
===========================

``vcffixup`` · *1 contributor · 1 version*

usage: vcffixup [file]
Count the allele frequencies across alleles
 present in each record in the VCF file. (Similar to vcftools --freq.)

Uses genotypes from the VCF file to correct AC (alternate allele count), AF (alternate allele frequency), NS (number of called), in the VCF records.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcffixup.versions import VcfFixUp_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcffixup_step",
           VcfFixUp_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcffixup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcffixup:

.. code-block:: bash

   # user inputs
   janis inputs vcffixup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run vcffixup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcffixup





Information
------------

:ID: ``vcffixup``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18


Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>  VCF output
======  ===========  ===============


Additional configuration (inputs)
---------------------------------

======  ============  ========  ==========  ===============
name    type          prefix      position  documentation
======  ============  ========  ==========  ===============
vcf     Gzipped<VCF>                     3
======  ============  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcffixup {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
     }
     command <<<
       set -e
       vcffixup \
         '~{vcf}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "shollizeck/vcflib:1.0.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'VcfLib: VcfFixUp'
   doc: |-
     usage: vcffixup [file]
     Count the allele frequencies across alleles
      present in each record in the VCF file. (Similar to vcftools --freq.)

     Uses genotypes from the VCF file to correct AC (alternate allele count), AF (alternate allele frequency), NS (number of called), in the VCF records.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/vcflib:1.0.1

   inputs:
   - id: vcf
     label: vcf
     type: File
     inputBinding:
       position: 3

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcffixup
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcffixup


