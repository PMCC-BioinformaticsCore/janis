:orphan:

Combine Variants
==================================

``combinevariants`` · *2 contributors · 1 version*


usage: combine_vcf.py [-h] -i I --columns COLUMNS -o O --type
                      {germline,somatic} [--regions REGIONS] [--normal NORMAL]
                      [--tumor TUMOR] [--priority PRIORITY [PRIORITY ...]]

Extracts and combines the information from germline / somatic vcfs into one

required arguments:
  -i I                  input vcfs, the priority of the vcfs will be based on
                        the order of the input. This parameter can be
                        specified more than once
  --columns COLUMNS     Columns to keep. This parameter can be specified more
                        than once
  -o O                  output vcf (unsorted)
  --type {germline,somatic}
                        must be either germline or somatic
  --regions REGIONS     Region file containing all the variants, used as
                        samtools mpileup
  --normal NORMAL       Sample id of germline vcf, or normal sample id of
                        somatic vcf
  --tumor TUMOR         tumor sample ID, required if inputs are somatic vcfs
  --priority PRIORITY [PRIORITY ...]
                        The priority of the callers, must match with the
                        callers in the source header

optional arguments:
  -h, --help            show this help message and exit



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.combinevariants.versions import CombineVariants_0_0_8

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "combinevariants_step",
           CombineVariants_0_0_8(
               vcfs=None,
               type=None,
           )
       )
       wf.output("out", source=combinevariants_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for combinevariants:

.. code-block:: bash

   # user inputs
   janis inputs combinevariants > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       type: <value>
       vcfs:
       - vcfs_0.vcf
       - vcfs_1.vcf




5. Run combinevariants with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       combinevariants





Information
------------

:ID: ``combinevariants``
:URL: `https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils <https://github.com/PMCC-BioinformaticsCore/scripts/tree/master/vcf_utils>`_
:Versions: 0.0.8
:Container: michaelfranklin/pmacutil:0.0.8
:Authors: Jiaan Yu, Michael Franklin
:Citations: None
:Created: 2019-03-25 00:00:00
:Updated: 2019-07-04 00:00:00


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  =======================  ==========  ==========  =============================================================================
name            type                     prefix      position    documentation
==============  =======================  ==========  ==========  =============================================================================
vcfs            Array<VCF>               -i                      input vcfs, the priority of the vcfs will be based on the order of the input
type            String                   --type                  germline | somatic
outputFilename  Optional<Filename>       -o
columns         Optional<Array<String>>  --columns               Columns to keep, seperated by space output vcf (unsorted)
normal          Optional<String>         --normal                Sample id of germline vcf, or normal sample id of somatic vcf
tumor           Optional<String>         --tumor                 tumor sample ID, required if inputs are somatic vcfs
priority        Optional<Integer>        --priority              The priority of the callers, must match with the callers in the source header
==============  =======================  ==========  ==========  =============================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task combinevariants {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? outputFilename
       Array[File] vcfs
       String type
       Array[String]? columns
       String? normal
       String? tumor
       Int? priority
     }
     command <<<
       set -e
       combine_vcf.py \
         -o '~{select_first([outputFilename, "generated.combined.vcf"])}' \
         ~{if length(vcfs) > 0 then "-i '" + sep("' -i '", vcfs) + "'" else ""} \
         --type '~{type}' \
         ~{if (defined(columns) && length(select_first([columns])) > 0) then "--columns '" + sep("','", select_first([columns])) + "'" else ""} \
         ~{if defined(normal) then ("--normal '" + normal + "'") else ""} \
         ~{if defined(tumor) then ("--tumor '" + tumor + "'") else ""} \
         ~{if defined(priority) then ("--priority " + priority) else ''}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/pmacutil:0.0.8"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.combined.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Combine Variants
   doc: |2

     usage: combine_vcf.py [-h] -i I --columns COLUMNS -o O --type
                           {germline,somatic} [--regions REGIONS] [--normal NORMAL]
                           [--tumor TUMOR] [--priority PRIORITY [PRIORITY ...]]

     Extracts and combines the information from germline / somatic vcfs into one

     required arguments:
       -i I                  input vcfs, the priority of the vcfs will be based on
                             the order of the input. This parameter can be
                             specified more than once
       --columns COLUMNS     Columns to keep. This parameter can be specified more
                             than once
       -o O                  output vcf (unsorted)
       --type {germline,somatic}
                             must be either germline or somatic
       --regions REGIONS     Region file containing all the variants, used as
                             samtools mpileup
       --normal NORMAL       Sample id of germline vcf, or normal sample id of
                             somatic vcf
       --tumor TUMOR         tumor sample ID, required if inputs are somatic vcfs
       --priority PRIORITY [PRIORITY ...]
                             The priority of the callers, must match with the
                             callers in the source header

     optional arguments:
       -h, --help            show this help message and exit

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/pmacutil:0.0.8

   inputs:
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.combined.vcf
     inputBinding:
       prefix: -o
   - id: vcfs
     label: vcfs
     doc: input vcfs, the priority of the vcfs will be based on the order of the input
     type:
       type: array
       inputBinding:
         prefix: -i
       items: File
     inputBinding: {}
   - id: type
     label: type
     doc: germline | somatic
     type: string
     inputBinding:
       prefix: --type
   - id: columns
     label: columns
     doc: Columns to keep, seperated by space output vcf (unsorted)
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --columns
       itemSeparator: ','
   - id: normal
     label: normal
     doc: Sample id of germline vcf, or normal sample id of somatic vcf
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --normal
   - id: tumor
     label: tumor
     doc: tumor sample ID, required if inputs are somatic vcfs
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --tumor
   - id: priority
     label: priority
     doc: The priority of the callers, must match with the callers in the source header
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --priority

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.combined.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: combine_vcf.py
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: combinevariants


