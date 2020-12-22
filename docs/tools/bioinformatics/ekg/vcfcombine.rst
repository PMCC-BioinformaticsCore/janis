:orphan:

VcfLib: VcfCombine
===============================

``vcfcombine`` · *1 contributor · 1 version*

usage: vcfcombine [vcf file] [vcf file] ...

 options:
-h --help	This text.
-r --region REGION	A region specifier of the form chrN:x-y to bound the merge

Combines VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfcombine.versions import VcfCombine_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfcombine_step",
           VcfCombine_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfcombine_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfcombine:

.. code-block:: bash

   # user inputs
   janis inputs vcfcombine > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf:
       - vcf_0.vcf
       - vcf_1.vcf




5. Run vcfcombine with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfcombine





Information
------------

:ID: ``vcfcombine``
:URL: `https://github.com/vcflib/vcflib <https://github.com/vcflib/vcflib>`_
:Versions: v1.0.1
:Container: shollizeck/vcflib:1.0.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-18
:Updated: 2019-10-18


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<File>  VCF output
======  ============  ===============


Additional configuration (inputs)
---------------------------------

======  ================  ========  ==========  ==========================================================
name    type              prefix      position  documentation
======  ================  ========  ==========  ==========================================================
vcf     Array<VCF>                           2
region  Optional<String>  -r                 1  A region specifier of the form chrN:x-y to bound the merge
======  ================  ========  ==========  ==========================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcfcombine {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] vcf
       String? region
     }
     command <<<
       set -e
       vcfcombine \
         ~{if defined(region) then ("-r '" + region + "'") else ""} \
         ~{if length(vcf) > 0 then "'" + sep("' '", vcf) + "'" else ""}
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
   label: 'VcfLib: VcfCombine'
   doc: |-
     usage: vcfcombine [vcf file] [vcf file] ...

      options:
     -h --help	This text.
     -r --region REGION	A region specifier of the form chrN:x-y to bound the merge

     Combines VCF files positionally, combining samples when sites and alleles are identical. Any number of VCF files may be combined. The INFO field and other columns are taken from one of the files which are combined when records in multiple files match. Alleles must have identical ordering to be combined into one record. If they do not, multiple records will be emitted

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/vcflib:1.0.1

   inputs:
   - id: vcf
     label: vcf
     type:
       type: array
       items: File
     inputBinding:
       position: 2
   - id: region
     label: region
     doc: A region specifier of the form chrN:x-y to bound the merge
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -r
       position: 1

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcfcombine
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcfcombine


