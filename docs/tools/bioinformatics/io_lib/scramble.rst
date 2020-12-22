:orphan:

scramble
========

``scramble`` · *1 contributor · 1 version*

scramble: streaming bam to cram compression


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.io_lib.scramble.versions import Scramble_1_14_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "scramble_step",
           Scramble_1_14_1_2(
               inputFilename=None,
               reference=None,
           )
       )
       wf.output("out", source=scramble_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for scramble:

.. code-block:: bash

   # user inputs
   janis inputs scramble > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputFilename: inputFilename.bam
       reference: reference.fasta




5. Run scramble with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       scramble





Information
------------

:ID: ``scramble``
:URL: `https://github.com/jkbonfield/io_lib/ <https://github.com/jkbonfield/io_lib/>`_
:Versions: 1.14.12
:Container: quay.io/biocontainers/staden_io_lib:1.14.12--h244ad75_0
:Authors: Matthias De Smet (@mattdsm)
:Citations: None
:Created: 2020-02-27
:Updated: 2020-02-27


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     stdout<CRAM>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

==========================  ==================  ========  ==========  =================================================
name                        type                prefix      position  documentation
==========================  ==================  ========  ==========  =================================================
inputFilename               BAM                                  200
reference                   FastaFai            -r                    Reference sequence file.
outputFilename              Optional<Filename>
range                       Optional<String>    -R                    Specifies the refseq:start-end range
maxBases                    Optional<Integer>   -b                    Max. bases per slice, default 5000000.
maxSequences                Optional<Integer>   -s                    Sequences per slice, default 10000.
maxSlicesPerContainer       Optional<Integer>   -S                    Slices per container, default 1.
embedReferenceSeuence       Optional<Boolean>   -e                    Embed reference sequence.
nonReferenceBaseEncoding    Optional<Boolean>   -x                    Non-reference based encoding.
multipleReferencesPerSlice  Optional<Boolean>   -M                    Use multiple references per slice.
generateTags                Optional<Boolean>   -m                    Generate MD and NM tags.
lzmaCompression             Optional<Boolean>   -Z                    Also compress using lzma
discardReadNames            Optional<Boolean>   -n                    Discard read names where possible.
preserveAuxTags             Optional<Boolean>   -P                    Preserve all aux tags (incl RG,NM,MD).
preserveAuxTagSizes         Optional<Boolean>   -p                    Preserve aux tag sizes ('i', 's', 'c').
noAddPG                     Optional<Boolean>   -q                    Don't add scramble @PG header line.
decodeStop                  Optional<Integer>   -N                    Stop decoding after 'integer' sequences.
threads                     Optional<Integer>   -t                    Number of threads. (default = 1)
enableQualityBinning        Optional<Integer>   -B                    Enable Illumina 8 quality-binning system (lossy).
==========================  ==================  ========  ==========  =================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task scramble {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File inputFilename
       File reference
       File reference_fai
       String? outputFilename
       String? range
       Int? maxBases
       Int? maxSequences
       Int? maxSlicesPerContainer
       Boolean? embedReferenceSeuence
       Boolean? nonReferenceBaseEncoding
       Boolean? multipleReferencesPerSlice
       Boolean? generateTags
       Boolean? lzmaCompression
       Boolean? discardReadNames
       Boolean? preserveAuxTags
       Boolean? preserveAuxTagSizes
       Boolean? noAddPG
       Int? decodeStop
       Int? threads
       Int? enableQualityBinning
     }
     command <<<
       set -e
       scramble \
         -r '~{reference}' \
         ~{if defined(range) then ("-R '" + range + "'") else ""} \
         ~{if defined(select_first([maxBases, 5000000])) then ("-b " + select_first([maxBases, 5000000])) else ''} \
         ~{if defined(select_first([maxSequences, 10000])) then ("-s " + select_first([maxSequences, 10000])) else ''} \
         ~{if defined(select_first([maxSlicesPerContainer, 1])) then ("-S " + select_first([maxSlicesPerContainer, 1])) else ''} \
         ~{if (defined(embedReferenceSeuence) && select_first([embedReferenceSeuence])) then "-e" else ""} \
         ~{if (defined(nonReferenceBaseEncoding) && select_first([nonReferenceBaseEncoding])) then "-x" else ""} \
         ~{if (defined(multipleReferencesPerSlice) && select_first([multipleReferencesPerSlice])) then "-M" else ""} \
         ~{if (defined(generateTags) && select_first([generateTags])) then "-m" else ""} \
         ~{if (defined(lzmaCompression) && select_first([lzmaCompression])) then "-Z" else ""} \
         ~{if (defined(discardReadNames) && select_first([discardReadNames])) then "-n" else ""} \
         ~{if (defined(preserveAuxTags) && select_first([preserveAuxTags])) then "-P" else ""} \
         ~{if (defined(preserveAuxTagSizes) && select_first([preserveAuxTagSizes])) then "-p" else ""} \
         ~{if (defined(noAddPG) && select_first([noAddPG])) then "-q" else ""} \
         ~{if defined(decodeStop) then ("-N " + decodeStop) else ''} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("-t " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(enableQualityBinning) then ("-B " + enableQualityBinning) else ''} \
         -I 'bam' \
         -O 'cram' \
         '-9' \
         -V '3.0' \
         '~{inputFilename}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/staden_io_lib:1.14.12--h244ad75_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
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
   label: scramble
   doc: 'scramble: streaming bam to cram compression'

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/staden_io_lib:1.14.12--h244ad75_0

   inputs:
   - id: inputFilename
     label: inputFilename
     type: File
     inputBinding:
       position: 200
   - id: reference
     label: reference
     doc: Reference sequence file.
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: -r
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.bam
   - id: range
     label: range
     doc: Specifies the refseq:start-end range
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -R
   - id: maxBases
     label: maxBases
     doc: Max. bases per slice, default 5000000.
     type: int
     default: 5000000
     inputBinding:
       prefix: -b
   - id: maxSequences
     label: maxSequences
     doc: Sequences per slice, default 10000.
     type: int
     default: 10000
     inputBinding:
       prefix: -s
   - id: maxSlicesPerContainer
     label: maxSlicesPerContainer
     doc: Slices per container, default 1.
     type: int
     default: 1
     inputBinding:
       prefix: -S
   - id: embedReferenceSeuence
     label: embedReferenceSeuence
     doc: Embed reference sequence.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -e
   - id: nonReferenceBaseEncoding
     label: nonReferenceBaseEncoding
     doc: Non-reference based encoding.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -x
   - id: multipleReferencesPerSlice
     label: multipleReferencesPerSlice
     doc: Use multiple references per slice.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -M
   - id: generateTags
     label: generateTags
     doc: Generate MD and NM tags.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -m
   - id: lzmaCompression
     label: lzmaCompression
     doc: Also compress using lzma
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -Z
   - id: discardReadNames
     label: discardReadNames
     doc: Discard read names where possible.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -n
   - id: preserveAuxTags
     label: preserveAuxTags
     doc: Preserve all aux tags (incl RG,NM,MD).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -P
   - id: preserveAuxTagSizes
     label: preserveAuxTagSizes
     doc: Preserve aux tag sizes ('i', 's', 'c').
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
   - id: noAddPG
     label: noAddPG
     doc: Don't add scramble @PG header line.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -q
   - id: decodeStop
     label: decodeStop
     doc: Stop decoding after 'integer' sequences.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -N
   - id: threads
     label: threads
     doc: Number of threads. (default = 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -t
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
   - id: enableQualityBinning
     label: enableQualityBinning
     doc: Enable Illumina 8 quality-binning system (lossy).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -B

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - scramble
   arguments:
   - prefix: -I
     position: 0
     valueFrom: bam
   - prefix: -O
     position: 0
     valueFrom: cram
   - position: 0
     valueFrom: '-9'
   - prefix: -V
     position: 0
     valueFrom: '3.0'

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: scramble


