:orphan:

GATK4: CreateSequenceDictionary
===============================================================

``Gatk4CreateSequenceDictionary`` · *1 contributor · 3 versions*

Creates a sequence dictionary for a reference sequence.  This tool creates a sequence dictionary file (with ".dict"
extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools.
The output file contains a header but no SAMRecords, and the header contains only sequence records.

The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

Usage example:

    java -jar picard.jar CreateSequenceDictionary \
        R=reference.fasta \
        O=reference.dict


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.createsequencedictionary.versions import Gatk4CreateSequenceDictionary_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4createsequencedictionary_step",
           Gatk4CreateSequenceDictionary_4_1_2(
               reference=None,
           )
       )
       wf.output("out", source=gatk4createsequencedictionary_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4CreateSequenceDictionary:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4CreateSequenceDictionary > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       reference: reference.fasta




5. Run Gatk4CreateSequenceDictionary with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4CreateSequenceDictionary





Information
------------

:ID: ``Gatk4CreateSequenceDictionary``
:URL: `https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard- <https://gatk.broadinstitute.org/hc/en-us/articles/360036509572-CreateSequenceDictionary-Picard->`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: TBD
:Created: 2020-02-14
:Updated: 2020-02-14


Outputs
-----------

======  ========  ======================================
name    type      documentation
======  ========  ======================================
out     FastDict  Output reference with ^.dict reference
======  ========  ======================================


Additional configuration (inputs)
---------------------------------

=================  =======================  ===========  ==========  ========================================================================================
name               type                     prefix       position    documentation
=================  =======================  ===========  ==========  ========================================================================================
reference          Fasta                    --REFERENCE              (-R) Input reference fasta or fasta.gz  Required.
javaOptions        Optional<Array<String>>
compression_level  Optional<Integer>                                 Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
=================  =======================  ===========  ==========  ========================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4CreateSequenceDictionary {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File reference
     }
     command <<<
       set -e
       cp -f '~{reference}' '.'
       gatk CreateSequenceDictionary \
         --java-options '-Xmx~{((select_first([runtime_memory, 2, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --REFERENCE '~{basename(reference)}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 2, 4])}G"
       preemptible: 2
     }
     output {
       File out = basename(reference)
       File out_dict = sub(sub(sub(basename(reference), "\\.fasta$", ".dict"), "\\.fna$", ".dict"), "\\.fa$", ".dict")
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: CreateSequenceDictionary'
   doc: |-
     Creates a sequence dictionary for a reference sequence.  This tool creates a sequence dictionary file (with ".dict"
     extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools.
     The output file contains a header but no SAMRecords, and the header contains only sequence records.

     The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

     Usage example:

         java -jar picard.jar CreateSequenceDictionary \
             R=reference.fasta \
             O=reference.dict

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.reference)
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.2.0

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: reference
     label: reference
     doc: (-R) Input reference fasta or fasta.gz  Required.
     type: File
     inputBinding:
       prefix: --REFERENCE

   outputs:
   - id: out
     label: out
     doc: Output reference with ^.dict reference
     type: File
     secondaryFiles:
     - pattern: ^.dict
     outputBinding:
       glob: $(inputs.reference.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - CreateSequenceDictionary
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 2, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4CreateSequenceDictionary


