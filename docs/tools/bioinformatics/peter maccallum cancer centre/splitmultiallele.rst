:orphan:

Split Multiple Alleles
=========================================

``SplitMultiAllele`` · *2 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.splitmultiallele import SplitMultiAllele

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "splitmultiallele_step",
           SplitMultiAllele(
               vcf=None,
               reference=None,
           )
       )
       wf.output("out", source=splitmultiallele_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SplitMultiAllele:

.. code-block:: bash

   # user inputs
   janis inputs SplitMultiAllele > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta
       vcf: vcf.vcf




5. Run SplitMultiAllele with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SplitMultiAllele





Information
------------

:ID: ``SplitMultiAllele``
:URL: *No URL to the documentation was provided*
:Versions: v0.5772
:Container: heuermh/vt
:Authors: Michael Franklin, Jiaan Yu
:Citations: None
:Created: 2019-01-18
:Updated: 2020-11-06


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
vcf             VCF                                    1
reference       FastaWithIndexes    -r                 4
outputFilename  Optional<Filename>  -o                 6
==============  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SplitMultiAllele {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
     }
     command <<<
       set -e
        \
         vt decompose -s \
         ~{vcf} \
         | vt normalize -n -q - \
         -r ~{reference} \
         -o ~{select_first([outputFilename, "~{basename(vcf, ".vcf")}.norm.vcf"])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "heuermh/vt"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(vcf, ".vcf")}.norm.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Split Multiple Alleles
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: heuermh/vt

   inputs:
   - id: vcf
     label: vcf
     type: File
     inputBinding:
       position: 1
       shellQuote: false
   - id: reference
     label: reference
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: -r
       position: 4
       shellQuote: false
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.norm.vcf
     inputBinding:
       prefix: -o
       position: 6
       valueFrom: $(inputs.vcf.basename.replace(/.vcf$/, "")).norm.vcf
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.vcf.basename.replace(/.vcf$/, "")).norm.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: 'vt decompose -s '
     shellQuote: false
   - position: 2
     valueFrom: '| vt normalize -n -q - '
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SplitMultiAllele


