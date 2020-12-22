:orphan:

Trim IUPAC Bases
============================

``trimIUPAC`` · *1 contributor · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.trimiupac.versions import TrimIUPAC_0_0_5

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "trimiupac_step",
           TrimIUPAC_0_0_5(
               vcf=None,
           )
       )
       wf.output("out", source=trimiupac_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for trimIUPAC:

.. code-block:: bash

   # user inputs
   janis inputs trimIUPAC > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run trimIUPAC with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       trimIUPAC





Information
------------

:ID: ``trimIUPAC``
:URL: *No URL to the documentation was provided*
:Versions: 0.0.5, 0.0.4
:Container: michaelfranklin/pmacutil:0.0.5
:Authors: Michael Franklin
:Citations: None
:Created: 2019-05-30
:Updated: 2019-12-08


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ======================================
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ======================================
vcf             VCF                                    0  The VCF to remove the IUPAC bases from
outputFilename  Optional<Filename>                     2
==============  ==================  ========  ==========  ======================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task trimIUPAC {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? outputFilename
     }
     command <<<
       set -e
       trimIUPAC.py \
         '~{vcf}' \
         '~{select_first([outputFilename, "generated.trimmed.vcf"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.0.5"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 1, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.trimmed.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Trim IUPAC Bases
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.5

   inputs:
   - id: vcf
     label: vcf
     doc: The VCF to remove the IUPAC bases from
     type: File
     inputBinding:
       position: 0
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.trimmed.vcf
     inputBinding:
       position: 2

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.trimmed.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: trimIUPAC.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: trimIUPAC


