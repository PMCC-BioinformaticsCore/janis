:orphan:

VcfLib: VcfStreamSort
=====================================

``vcfstreamsort`` · *1 contributor · 1 version*

usage: vcfallelicprimitives [options] [file]

options:
	-m, --use-mnps	Retain MNPs as separate events (default: false)
	-t, --tag-parsed FLAG	Tag records which are split apart of a complex allele with this flag


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfstreamsort.versions import VcfStreamSort_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfstreamsort_step",
           VcfStreamSort_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfstreamsort_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfstreamsort:

.. code-block:: bash

   # user inputs
   janis inputs vcfstreamsort > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run vcfstreamsort with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfstreamsort





Information
------------

:ID: ``vcfstreamsort``
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

============  =================  ========  ==========  =======================================
name          type               prefix      position  documentation
============  =================  ========  ==========  =======================================
vcf           VCF                                   3
inMemoryFlag  Optional<Boolean>  -a                    load all sites and then sort in memory
windowSize    Optional<Integer>  -w                    number of sites to sort (default 10000)
============  =================  ========  ==========  =======================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcfstreamsort {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       Boolean? inMemoryFlag
       Int? windowSize
     }
     command <<<
       set -e
       vcfstreamsort \
         ~{if select_first([inMemoryFlag, false]) then "-a" else ""} \
         ~{if defined(windowSize) then ("-w " + windowSize) else ''} \
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
   label: 'VcfLib: VcfStreamSort'
   doc: |-
     usage: vcfallelicprimitives [options] [file]

     options:
     	-m, --use-mnps	Retain MNPs as separate events (default: false)
     	-t, --tag-parsed FLAG	Tag records which are split apart of a complex allele with this flag

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
   - id: inMemoryFlag
     label: inMemoryFlag
     doc: load all sites and then sort in memory
     type: boolean
     default: false
     inputBinding:
       prefix: -a
   - id: windowSize
     label: windowSize
     doc: number of sites to sort (default 10000)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -w

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcfstreamsort
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcfstreamsort


