:orphan:

VcfLib: VcfAllelicPrimitives
===================================================

``vcfallelicprimitives`` · *1 contributor · 1 version*

usage: vcfallelicprimitives [options] [file]

options:
	-m, --use-mnps	Retain MNPs as separate events (default: false)
	-t, --tag-parsed FLAG	Tag records which are split apart of a complex allele with this flag


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcflib.vcfallelicprimitives.versions import VcfAllelicPrimitives_1_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcfallelicprimitives_step",
           VcfAllelicPrimitives_1_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=vcfallelicprimitives_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vcfallelicprimitives:

.. code-block:: bash

   # user inputs
   janis inputs vcfallelicprimitives > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run vcfallelicprimitives with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vcfallelicprimitives





Information
------------

:ID: ``vcfallelicprimitives``
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

============  =================  ========  ==========  =======================================================================================================================================================================================================================================
name          type               prefix      position  documentation
============  =================  ========  ==========  =======================================================================================================================================================================================================================================
vcf           Gzipped<VCF>                          3
useMnpsFlag   Optional<Boolean>  -m                    Retain MNPs as separate events (default: false)
tagParsed     Optional<String>   -t                    Tag records which are split apart of a complex allele with this flag
keepInfoFlag  Optional<Boolean>  -k                    Maintain site and allele-level annotations when decomposing. Note that in many cases, such as multisample VCFs, these won't be valid post-decomposition.  For biallelic loci in single-sample VCFs, they should be usable with caution.
keepGenoFlag  Optional<Boolean>  -g                    Maintain genotype-level annotations when decomposing.  Similar caution should be used for this as for --keep-info.
maxLength     Optional<Integer>  -L                    Do not manipulate records in which either the ALT or REF is longer than LEN (default: 200).
============  =================  ========  ==========  =======================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vcfallelicprimitives {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       Boolean? useMnpsFlag
       String? tagParsed
       Boolean? keepInfoFlag
       Boolean? keepGenoFlag
       Int? maxLength
     }
     command <<<
       set -e
       vcfallelicprimitives \
         ~{if select_first([useMnpsFlag, false]) then "-m" else ""} \
         ~{if defined(tagParsed) then ("-t '" + tagParsed + "'") else ""} \
         ~{if (defined(keepInfoFlag) && select_first([keepInfoFlag])) then "-k" else ""} \
         ~{if (defined(keepGenoFlag) && select_first([keepGenoFlag])) then "-g" else ""} \
         ~{if defined(maxLength) then ("-L " + maxLength) else ''} \
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
   label: 'VcfLib: VcfAllelicPrimitives'
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
   - id: useMnpsFlag
     label: useMnpsFlag
     doc: 'Retain MNPs as separate events (default: false)'
     type: boolean
     default: false
     inputBinding:
       prefix: -m
   - id: tagParsed
     label: tagParsed
     doc: Tag records which are split apart of a complex allele with this flag
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -t
   - id: keepInfoFlag
     label: keepInfoFlag
     doc: |-
       Maintain site and allele-level annotations when decomposing. Note that in many cases, such as multisample VCFs, these won't be valid post-decomposition.  For biallelic loci in single-sample VCFs, they should be usable with caution.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -k
   - id: keepGenoFlag
     label: keepGenoFlag
     doc: |-
       Maintain genotype-level annotations when decomposing.  Similar caution should be used for this as for --keep-info.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -g
   - id: maxLength
     label: maxLength
     doc: |-
       Do not manipulate records in which either the ALT or REF is longer than LEN (default: 200).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -L

   outputs:
   - id: out
     label: out
     doc: VCF output
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand: vcfallelicprimitives
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vcfallelicprimitives


