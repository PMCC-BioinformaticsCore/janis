:orphan:

VcfLib: Vcf Length
==============================

``vcflength`` · *1 contributor · 1 version*

Adds the length of the variant record (in [-/+]) relative to the reference allele to each VCF record.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcflength.versions import VcfLength_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcflength_step",
           VcfLength_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcflength_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcflength:

.. code-block:: bash

   # user inputs
   janis inputs vcflength > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcflength with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcflength





Information
------------

:ID: ``vcflength``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Michael Franklin
:Citations: None
:Created: 2020-06-04
:Updated: 2020-06-04


Outputs
-----------

======  ===========  ==========================================================
name    type         documentation
======  ===========  ==========================================================
out     stdout<VCF>  VCF with length of the variant record added to each record
======  ===========  ==========================================================


Additional configuration (inputs)
---------------------------------

======  ======  ========  ==========  ========================================================================
name    type    prefix      position  documentation
======  ======  ========  ==========  ========================================================================
vcf     VCF                        1  VCF to add length of variant record relative to the reference allele to.
======  ======  ========  ==========  ========================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcflength {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
     }
     command <<<
       set -e
       vcflength \
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
   label: 'VcfLib: Vcf Length'
   doc: |-
     Adds the length of the variant record (in [-/+]) relative to the reference allele to each VCF record.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/vcflib:1.0.1

   inputs:
   - id: vcf
     label: vcf
     doc: VCF to add length of variant record relative to the reference allele to.
     type: File
     inputBinding:
       position: 1

   outputs:
   - id: out
     label: out
     doc: VCF with length of the variant record added to each record
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcflength
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcflength


