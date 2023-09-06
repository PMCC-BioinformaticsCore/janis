:orphan:

TransposeTsv
============

``TransposeTsv`` · *0 contributors · 1 version*

transpose tsv file


Quickstart
-----------

    .. code-block:: python

       from janis_unix.tools.transposetsv import TransposeTsv

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "transposetsv_step",
           TransposeTsv(
               inp_tsv=None,
           )
       )
       wf.output("out", source=transposetsv_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

4. Generate user input files for TransposeTsv:

.. code-block:: bash

   # user inputs
   janis inputs TransposeTsv > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp_tsv: inp_tsv.tsv




5. Run TransposeTsv with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       TransposeTsv

.. note::

   You can use `janis prepare <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ to improve setting up your files for this CommandTool. See `this guide <https://janis.readthedocs.io/en/latest/references/prepare.html>`_ for more information about Janis Prepare.

   .. code-block:: text

      OUTPUT_DIR="<output-dir>"
      janis prepare \
          --inputs inputs.yaml \
          --output-dir $OUTPUT_DIR \
          TransposeTsv

      # Run script that Janis automatically generates
      sh $OUTPUT_DIR/run.sh











Information
------------

:ID: ``TransposeTsv``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Container: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
inp_tsv         tsv                                    2
outputFilename  Optional<Filename>                     4
==============  ==================  ========  ==========  ===============

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task TransposeTsv {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disk
       File inp_tsv
       String? outputFilename
     }

     command <<<
       set -e
        \
         echo 'BEGIN { FS=OFS="	" }
   { 
       for (rowNr=1;rowNr<=NF;rowNr++) {
           cell[rowNr,NR] = $rowNr
       }
       maxRows = (NF > maxRows ? NF : maxRows)
       maxCols = NR
   }
   END {
       for (rowNr=1;rowNr<=maxRows;rowNr++) {
           for (colNr=1;colNr<=maxCols;colNr++) {
               printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
           }
       }
   }' > tst.awk; \
         awk -f tst.awk \
         '~{inp_tsv}' \
         > \
         '~{select_first([outputFilename, "~{basename(inp_tsv, ".tsv")}.transposed.tsv"])}'
     >>>

     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disk, 20])} SSD"
       docker: "ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }

     output {
       File out = select_first([outputFilename, "~{basename(inp_tsv, ".tsv")}.transposed.tsv"])
     }

   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: TransposeTsv

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: ubuntu@sha256:1d7b639619bdca2d008eca2d5293e3c43ff84cbee597ff76de3b7a7de3e84956

   inputs:
   - id: inp_tsv
     label: inp_tsv
     type: File
     inputBinding:
       position: 2
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.transposed.tsv
     inputBinding:
       position: 4
       valueFrom: $(inputs.inp_tsv.basename.replace(/.tsv$/, "")).transposed.tsv

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.inp_tsv.basename.replace(/.tsv$/, "")).transposed.tsv
       loadContents: false
   stdout: _stdout
   stderr: _stderr
   arguments:
   - position: 0
     valueFrom: |-
       echo 'BEGIN { FS=OFS="	" }
       { 
           for (rowNr=1;rowNr<=NF;rowNr++) {
               cell[rowNr,NR] = $rowNr
           }
           maxRows = (NF > maxRows ? NF : maxRows)
           maxCols = NR
       }
       END {
           for (rowNr=1;rowNr<=maxRows;rowNr++) {
               for (colNr=1;colNr<=maxCols;colNr++) {
                   printf "%s%s", cell[rowNr,colNr], (colNr < maxCols ? OFS : ORS)
               }
           }
       }' > tst.awk;
     shellQuote: false
   - position: 1
     valueFrom: 'awk -f tst.awk '
     shellQuote: false
   - position: 3
     valueFrom: '>'
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: TransposeTsv


