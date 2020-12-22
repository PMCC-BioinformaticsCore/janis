:orphan:

VcfLib: Vcf Random Sampling
=============================================

``vcfrandomsample`` · *1 contributor · 1 version*

usage: vcfrandomsample [options] [<vcf file>]

options:
	-r, --rate RATE 	base sampling probability per locus
	-s, --scale-by KEY\scale sampling likelihood by this Float info field
	-p, --random-seed N	use this random seed


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfrandomsample.versions import VcfRandomSample_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfrandomsample_step",
           VcfRandomSample_1_0_1(
               vcf=None,
               rate=None,
               seed=None,
           )
       )
       wf.output("out", source=vcfrandomsample_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfrandomsample:

.. code-block:: bash

   # user inputs
   janis inputs vcfrandomsample > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       rate: 0.0
       seed: 0
       vcf: vcf.vcf.gz




5. Run vcfrandomsample with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfrandomsample





Information
------------

:ID: ``vcfrandomsample``
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

=======  ================  ========  ==========  ==================================================
name     type              prefix      position  documentation
=======  ================  ========  ==========  ==================================================
vcf      Gzipped<VCF>                         3
rate     Float             -t                    base sampling probability per locus
seed     Integer           -p                    use this random seed
scaleBy  Optional<String>  -s                    scale sampling likelihood by this Float info field
=======  ================  ========  ==========  ==================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcfrandomsample {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       Float rate
       String? scaleBy
       Int seed
     }
     command <<<
       set -e
       vcfrandomsample \
         -t ~{rate} \
         ~{if defined(scaleBy) then ("-s '" + scaleBy + "'") else ""} \
         -p ~{seed} \
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
   label: 'VcfLib: Vcf Random Sampling'
   doc: |-
     usage: vcfrandomsample [options] [<vcf file>]

     options:
     	-r, --rate RATE 	base sampling probability per locus
     	-s, --scale-by KEY\scale sampling likelihood by this Float info field
     	-p, --random-seed N	use this random seed

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
   - id: rate
     label: rate
     doc: base sampling probability per locus
     type: float
     inputBinding:
       prefix: -t
   - id: scaleBy
     label: scaleBy
     doc: scale sampling likelihood by this Float info field
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -s
   - id: seed
     label: seed
     doc: use this random seed
     type: int
     inputBinding:
       prefix: -p

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcfrandomsample
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcfrandomsample


