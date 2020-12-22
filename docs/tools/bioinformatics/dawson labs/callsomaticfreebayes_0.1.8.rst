:orphan:

Call Somatic Variants from freebayes
===========================================================

``callSomaticFreeBayes`` · *1 contributor · 1 version*

Usage: callSomaticFreeBayes.R [options]



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.callsomaticfreebayes.latest import CallSomaticFreeBayesLatest

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "callsomaticfreebayes_step",
           CallSomaticFreeBayesLatest(
               vcf=None,
           )
       )
       wf.output("out", source=callsomaticfreebayes_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for callSomaticFreeBayes:

.. code-block:: bash

   # user inputs
   janis inputs callSomaticFreeBayes > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run callSomaticFreeBayes with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       callSomaticFreeBayes





Information
------------

:ID: ``callSomaticFreeBayes``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.8
:Container: shollizeck/dawsontoolkit:0.1.8.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-19
:Updated: 2019-10-25


Outputs
-----------

======  ======  =================
name    type    documentation
======  ======  =================
out     VCF     To determine type
======  ======  =================


Additional configuration (inputs)
---------------------------------

================  ==================  ========  ==========  ================================================================
name              type                prefix    position    documentation
================  ==================  ========  ==========  ================================================================
vcf               VCF                 -i                    input vcf
normalSampleName  Optional<String>    -n                    the normal sample name in the vcf (default: first sample in vcf)
outputFilename    Optional<Filename>  -o                    output file name (default: STDOUT)
================  ==================  ========  ==========  ================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task callSomaticFreeBayes {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? normalSampleName
       String? outputFilename
     }
     command <<<
       set -e
       callSomaticFreeBayes.R \
         -i '~{vcf}' \
         ~{if defined(normalSampleName) then ("-n '" + normalSampleName + "'") else ""} \
         -o '~{select_first([outputFilename, "generated.vcf"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "shollizeck/dawsontoolkit:0.1.8.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 12, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Call Somatic Variants from freebayes
   doc: "Usage: callSomaticFreeBayes.R [options]\n"

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: shollizeck/dawsontoolkit:0.1.8.1

   inputs:
   - id: vcf
     label: vcf
     doc: input vcf
     type: File
     inputBinding:
       prefix: -i
   - id: normalSampleName
     label: normalSampleName
     doc: 'the normal sample name in the vcf (default: first sample in vcf)'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -n
   - id: outputFilename
     label: outputFilename
     doc: 'output file name (default: STDOUT)'
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: -o

   outputs:
   - id: out
     label: out
     doc: To determine type
     type: File
     outputBinding:
       glob: generated.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: callSomaticFreeBayes.R
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: callSomaticFreeBayes


