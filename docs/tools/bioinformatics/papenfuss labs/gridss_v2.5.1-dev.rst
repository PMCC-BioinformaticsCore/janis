:orphan:

Gridss
===============

``gridss`` · *1 contributor · 5 versions*

GRIDSS: the Genomic Rearrangement IDentification Software Suite

GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
break-end assembly, split read, and read pair evidence.

GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
any step (such as split read identification) can be replaced by another implementation that also outputs
using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
pipeline designed for interoperability with other tools.

If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.papenfuss.gridss.gridss import Gridss_2_5_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gridss_step",
           Gridss_2_5_1(
               bams=None,
               reference=None,
           )
       )
       wf.output("out", source=gridss_step.out)
       wf.output("assembly", source=gridss_step.assembly)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for gridss:

.. code-block:: bash

   # user inputs
   janis inputs gridss > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       reference: reference.fasta




5. Run gridss with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       gridss





Information
------------

:ID: ``gridss``
:URL: `https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation <https://github.com/PapenfussLab/gridss/wiki/GRIDSS-Documentation>`_
:Versions: v2.9.4, v2.6.2, v2.5.1-dev, v2.4.0, v2.2.3
:Container: michaelfranklin/gridss:2.5.1-dev2
:Authors: Michael Franklin
:Citations: Daniel L. Cameron, Jan Schröder, Jocelyn Sietsma Penington, Hongdo Do, Ramyar Molania, Alexander Dobrovic, Terence P. Speed and Anthony T. Papenfuss. GRIDSS: sensitive and specific genomic rearrangement detection using positional de Bruijn graph assembly. Genome Research, 2017 doi: 10.1101/gr.222109.117
:DOI: 10.1101/gr.222109.117
:Created: 2019-06-19
:Updated: 2019-08-20


Outputs
-----------

========  ======  ===============
name      type    documentation
========  ======  ===============
out       VCF
assembly  BAM
========  ======  ===============


Additional configuration (inputs)
---------------------------------

================  ==================  ============  ==========  ===============
name              type                prefix          position  documentation
================  ==================  ============  ==========  ===============
bams              Array<IndexedBam>                         10
reference         FastaWithIndexes    --reference            1
outputFilename    Optional<Filename>  --output               2
assemblyFilename  Optional<Filename>  --assembly             3
threads           Optional<Integer>   --threads
blacklist         Optional<bed>       --blacklist            4
tmpdir            Optional<String>    --workingdir
================  ==================  ============  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task gridss {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] bams
       Array[File] bams_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       String? assemblyFilename
       Int? threads
       File? blacklist
       String? tmpdir
     }
     command <<<
       set -e
       gridss.sh \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("--threads " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(select_first([tmpdir, "./TMP"])) then ("--workingdir '" + select_first([tmpdir, "./TMP"]) + "'") else ""} \
         --reference '~{reference}' \
         --output '~{select_first([outputFilename, "generated.svs.vcf"])}' \
         --assembly '~{select_first([assemblyFilename, "generated.assembled.bam"])}' \
         ~{if defined(blacklist) then ("--blacklist '" + blacklist + "'") else ""} \
         ~{if length(bams) > 0 then "'" + sep("' '", bams) + "'" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 8, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/gridss:2.5.1-dev2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 31, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.svs.vcf"])
       File assembly = select_first([assemblyFilename, "generated.assembled.bam"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Gridss
   doc: |
     GRIDSS: the Genomic Rearrangement IDentification Software Suite

     GRIDSS is a module software suite containing tools useful for the detection of genomic rearrangements.
     GRIDSS includes a genome-wide break-end assembler, as well as a structural variation caller for Illumina
     sequencing data. GRIDSS calls variants based on alignment-guided positional de Bruijn graph genome-wide
     break-end assembly, split read, and read pair evidence.

     GRIDSS makes extensive use of the standard tags defined by SAM specifications. Due to the modular design,
     any step (such as split read identification) can be replaced by another implementation that also outputs
     using the standard tags. It is hoped that GRIDSS can serve as an exemplar modular structural variant
     pipeline designed for interoperability with other tools.

     If you have any trouble running GRIDSS, please raise an issue using the Issues tab above. Based on feedback
     from users, a user guide will be produced outlining common workflows, pitfalls, and use cases.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/gridss:2.5.1-dev2

   inputs:
   - id: bams
     label: bams
     type:
       type: array
       items: File
     inputBinding:
       position: 10
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
       prefix: --reference
       position: 1
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.svs.vcf
     inputBinding:
       prefix: --output
       position: 2
   - id: assemblyFilename
     label: assemblyFilename
     type:
     - string
     - 'null'
     default: generated.assembled.bam
     inputBinding:
       prefix: --assembly
       position: 3
   - id: threads
     label: threads
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
       valueFrom: |-
         $([inputs.runtime_cpu, 8, 1].filter(function (inner) { return inner != null })[0])
   - id: blacklist
     label: blacklist
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --blacklist
       position: 4
   - id: tmpdir
     label: tmpdir
     type: string
     default: ./TMP
     inputBinding:
       prefix: --workingdir

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.svs.vcf
       loadContents: false
   - id: assembly
     label: assembly
     type: File
     outputBinding:
       glob: generated.assembled.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: gridss.sh
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: gridss


