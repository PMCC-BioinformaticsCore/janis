:orphan:

FixUp FreeBayes MNPs
=========================================

``FixUpFreeBayesMNPs`` · *1 contributor · 1 version*

Usage: fixupFreeBayesMNPs.R [options]



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.fixupfreebayesmnps.versions import FixUpFreeBayesMNPs_0_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "fixupfreebayesmnps_step",
           FixUpFreeBayesMNPs_0_1(
               vcf=None,
           )
       )
       wf.output("out", source=fixupfreebayesmnps_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for FixUpFreeBayesMNPs:

.. code-block:: bash

   # user inputs
   janis inputs FixUpFreeBayesMNPs > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run FixUpFreeBayesMNPs with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FixUpFreeBayesMNPs

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          FixUpFreeBayesMNPs

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``FixUpFreeBayesMNPs``
:URL: *No URL to the documentation was provided*
:Versions: 0.2
:Container: shollizeck/dawsontoolkit:0.2
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2021-06-03
:Updated: 2021-06-03


Outputs
-----------

======  ============  =================
name    type          documentation
======  ============  =================
out     Gzipped<VCF>  To determine type
======  ============  =================


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============================================
name            type                prefix    position    documentation
==============  ==================  ========  ==========  ===============================================
vcf             VCF                 -i                    input vcf
outputFilename  Optional<Filename>  -o                    output file name (default: reassembled.vcf.bgz)
uncompressed    Optional<Boolean>   -o                    output file name (default: reassembled.vcf.bgz)
==============  ==================  ========  ==========  ===============================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task FixUpFreeBayesMNPs {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File vcf
       String? outputFilename
       Boolean? uncompressed
     }

     command <<<
       set -e
       fixupFreeBayesMNPs.R \
         -i '~{vcf}' \
         -o '~{select_first([outputFilename, "generated.vcf"])}' \
         ~{if (defined(uncompressed) && select_first([uncompressed])) then "-o" else ""}
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "shollizeck/dawsontoolkit:0.2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 12, 4])}G"
       preemptible: 2
     }

     output {
       File out = (select_first([outputFilename, "generated.vcf"]) + ".bgz")
       File out_tbi = (select_first([outputFilename, "generated.vcf"]) + ".bgz") + ".tbi"
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: FixUp FreeBayes MNPs

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/dawsontoolkit:0.2

   inputs:
   - id: vcf
     label: vcf
     doc: input vcf
     type: File
     inputBinding:
       prefix: -i
   - id: outputFilename
     label: outputFilename
     doc: 'output file name (default: reassembled.vcf.bgz)'
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: -o
   - id: uncompressed
     label: uncompressed
     doc: 'output file name (default: reassembled.vcf.bgz)'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -o

   outputs:
   - id: out
     label: out
     doc: To determine type
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $((inputs.outputFilename + ".bgz"))
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: fixupFreeBayesMNPs.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: FixUpFreeBayesMNPs


