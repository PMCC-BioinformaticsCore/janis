:orphan:

Split Multiple Alleles and Normalise Vcf
=======================================================================

``SplitMultiAlleleNormaliseVcf`` · *1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.splitmultiallele_normalistvcf import SplitMultiAlleleNormaliseVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "splitmultiallelenormalisevcf_step",
           SplitMultiAlleleNormaliseVcf(
               reference=None,
           )
       )
       wf.output("out", source=splitmultiallelenormalisevcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SplitMultiAlleleNormaliseVcf:

.. code-block:: bash

   # user inputs
   janis inputs SplitMultiAlleleNormaliseVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run SplitMultiAlleleNormaliseVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SplitMultiAlleleNormaliseVcf





Information
------------

:ID: ``SplitMultiAlleleNormaliseVcf``
:URL: *No URL to the documentation was provided*
:Versions: v0.5772
:Container: heuermh/vt
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-06-04
:Updated: 2020-06-04


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==================  ======================  ========  ==========  ===============
name                type                    prefix      position  documentation
==================  ======================  ========  ==========  ===============
reference           FastaWithIndexes        -r                 4
vcf                 Optional<VCF>                              1
compressedTabixVcf  Optional<Gzipped<VCF>>                     1
compressedVcf       Optional<Gzipped<VCF>>                     1
outputFilename      Optional<Filename>      -o                 6
==================  ======================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SplitMultiAlleleNormaliseVcf {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File? vcf
       File? compressedTabixVcf
       File? compressedTabixVcf_tbi
       File? compressedVcf
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
         ~{compressedTabixVcf} \
         ~{compressedVcf} \
         | vt normalize -n -q - \
         -r ~{reference} \
         -o ~{select_first([outputFilename, "generated.norm.vcf"])}
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
       File out = select_first([outputFilename, "generated.norm.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Split Multiple Alleles and Normalise Vcf
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: heuermh/vt

   inputs:
   - id: vcf
     label: vcf
     type:
     - File
     - 'null'
     inputBinding:
       position: 1
       shellQuote: false
   - id: compressedTabixVcf
     label: compressedTabixVcf
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       position: 1
       shellQuote: false
   - id: compressedVcf
     label: compressedVcf
     type:
     - File
     - 'null'
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
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.norm.vcf
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
   id: SplitMultiAlleleNormaliseVcf


