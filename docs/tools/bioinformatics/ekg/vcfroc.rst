:orphan:

VcfLib: Vcf ROC generator
==================================

``vcfroc`` · *1 contributor · 1 version*

usage: vcfroc [options] [<vcf file>]

options:
	-t, --truth-vcf FILE	use this VCF as ground truth for ROC generation
	-w, --window-size N       compare records up to this many bp away (default 30)
	-r, --reference FILE	FASTA reference file


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfroc.versions import VcfRoc_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfroc_step",
           VcfRoc_1_0_1(
               vcf=None,
               truth=None,
               reference=None,
           )
       )
       wf.output("out", source=vcfroc_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfroc:

.. code-block:: bash

   # user inputs
   janis inputs vcfroc > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta
       truth: truth.vcf.gz
       vcf: vcf.vcf.gz




5. Run vcfroc with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfroc





Information
------------

:ID: ``vcfroc``
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

==========  =================  ========  ==========  ====================================================
name        type               prefix      position  documentation
==========  =================  ========  ==========  ====================================================
vcf         Gzipped<VCF>                          3
truth       Gzipped<VCF>       -t                    use this VCF as ground truth for ROC generation
reference   FastaWithIndexes   -r                    FASTA reference file
windowSize  Optional<Integer>  -w                    compare records up to this many bp away (default 30)
==========  =================  ========  ==========  ====================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcfroc {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       File truth
       Int? windowSize
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
     }
     command <<<
       set -e
       vcfroc \
         -t '~{truth}' \
         ~{if defined(select_first([windowSize, 30])) then ("-w " + select_first([windowSize, 30])) else ''} \
         -r '~{reference}' \
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
   label: 'VcfLib: Vcf ROC generator'
   doc: |-
     usage: vcfroc [options] [<vcf file>]

     options:
     	-t, --truth-vcf FILE	use this VCF as ground truth for ROC generation
     	-w, --window-size N       compare records up to this many bp away (default 30)
     	-r, --reference FILE	FASTA reference file

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
   - id: truth
     label: truth
     doc: use this VCF as ground truth for ROC generation
     type: File
     inputBinding:
       prefix: -t
   - id: windowSize
     label: windowSize
     doc: compare records up to this many bp away (default 30)
     type: int
     default: 30
     inputBinding:
       prefix: -w
   - id: reference
     label: reference
     doc: FASTA reference file
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

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcfroc
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcfroc


