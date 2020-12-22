:orphan:

Sequenza: seqz binning
====================================

``SeqzBinning`` · *2 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.sequenza.seqz_binning.versions import SequenzaBinning_3_0_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "seqzbinning_step",
           SequenzaBinning_3_0_0(
               seqz=None,
               window=None,
           )
       )
       wf.output("out", source=seqzbinning_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SeqzBinning:

.. code-block:: bash

   # user inputs
   janis inputs SeqzBinning > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       seqz: seqz
       window: 0




5. Run SeqzBinning with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SeqzBinning





Information
------------

:ID: ``SeqzBinning``
:URL: `http://www.cbs.dtu.dk/biotools/sequenza/ <http://www.cbs.dtu.dk/biotools/sequenza/>`_
:Versions: 3.0.0, 2.2.0.9000
:Container: sequenza/sequenza:3.0.0
:Authors: mumbler, evanwehi
:Citations: None
:Created: 2019-12-16
:Updated: 2019-12-16


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============


Additional configuration (inputs)
---------------------------------

===============  ==================  ========  ==========  ===================================================================
name             type                prefix      position  documentation
===============  ==================  ========  ==========  ===================================================================
seqz             File                --seqz             2  A seqz file.
window           Integer             --window           4  Window size used for binning the original seqz file. Default is 50.
output_filename  Optional<Filename>  -o                 6  Output file "-" for STDOUT
===============  ==================  ========  ==========  ===================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SeqzBinning {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File seqz
       Int window
       String? output_filename
     }
     command <<<
       set -e
       sequenza-utils seqz_binning \
         --seqz '~{seqz}' \
         --window ~{window} \
         -o '~{select_first([output_filename, "generated.gz"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "sequenza/sequenza:3.0.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([output_filename, "generated.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'Sequenza: seqz binning'
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: sequenza/sequenza:3.0.0

   inputs:
   - id: seqz
     label: seqz
     doc: A seqz file.
     type: File
     inputBinding:
       prefix: --seqz
       position: 2
   - id: window
     label: window
     doc: Window size used for binning the original seqz file. Default is 50.
     type: int
     inputBinding:
       prefix: --window
       position: 4
   - id: output_filename
     label: output_filename
     doc: Output file "-" for STDOUT
     type:
     - string
     - 'null'
     default: generated.gz
     inputBinding:
       prefix: -o
       position: 6

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - sequenza-utils
   - seqz_binning
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SeqzBinning


