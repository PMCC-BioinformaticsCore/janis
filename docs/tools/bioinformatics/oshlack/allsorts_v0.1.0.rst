:orphan:

Allsorts
========

``Allsorts`` · *2 contributors · 1 version*

usage: ALLSorts [-h] -samples SAMPLES [-labels LABELS]
                [-destination DESTINATION] [-test] [-train]
                [-model_dir MODEL_DIR] [-njobs NJOBS] [-cv CV] [-verbose]
                [-comparison] [-force] [-parents]
ALLSorts CLI



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.oshlack.allsorts.versions import AllSorts_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "allsorts_step",
           AllSorts_0_1_0(
               samples=None,
           )
       )
       wf.output("out_predictions", source=allsorts_step.out_predictions)
       wf.output("out_probabilities", source=allsorts_step.out_probabilities)
       wf.output("out_distributions", source=allsorts_step.out_distributions)
       wf.output("out_waterfalls", source=allsorts_step.out_waterfalls)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Allsorts:

.. code-block:: bash

   # user inputs
   janis inputs Allsorts > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       samples: samples.csv




5. Run Allsorts with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Allsorts





Information
------------

:ID: ``Allsorts``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Container: breons/allsorts:0.1.0
:Authors: Michael Franklin, Jiaan Yu
:Citations: None
:Created: 2020-09-02
:Updated: 2020-09-11


Outputs
-----------

=================  ======  ===============
name               type    documentation
=================  ======  ===============
out_predictions    csv
out_probabilities  csv
out_distributions  File
out_waterfalls     File
=================  ======  ===============


Additional configuration (inputs)
---------------------------------

===========  =================  ============  ==========  =============================================================================================================================================================================================================================================================
name         type               prefix        position    documentation
===========  =================  ============  ==========  =============================================================================================================================================================================================================================================================
samples      csv                -samples                  (-s)  Path to samples (rows) x genes (columns) csv file representing a raw counts matrix. Note: hg19 only supported currently, use other references at own risk.
labels       Optional<csv>      -labels                   (-l)  (Optional) Path to samples true labels. CSV with samples (rows) x [sample id, label] (cols). This will enable re-labelling mode. Note: labels must reflect naming conventions used within this tool. View the ALLSorts GitHub Wiki for further details.
destination  Optional<String>   -destination              (-d)  Path to where you want the final report to be saved.
verbose      Optional<Boolean>  -verbose                  (-v) (flag, default=False) Verbose. Print stage progress.
comparison   Optional<Boolean>  -comparison               Rebuild comparisons for labelled visualisations.
force        Optional<Boolean>  -force                    (-f) (flag, default=False) Force. Bypass warnings without user confirmation.
parents      Optional<Boolean>  -parents                  (-p) Include parent meta-subtypes in predictions. Note: This may remove previously unclassified samples.
===========  =================  ============  ==========  =============================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Allsorts {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File samples
       File? labels
       String? destination
       Boolean? verbose
       Boolean? comparison
       Boolean? force
       Boolean? parents
     }
     command <<<
       set -e
       ALLSorts \
         -samples '~{samples}' \
         ~{if defined(labels) then ("-labels '" + labels + "'") else ""} \
         ~{if defined(select_first([destination, "."])) then ("-destination '" + select_first([destination, "."]) + "'") else ""} \
         ~{if (defined(verbose) && select_first([verbose])) then "-verbose" else ""} \
         ~{if (defined(comparison) && select_first([comparison])) then "-comparison" else ""} \
         ~{if (defined(force) && select_first([force])) then "-force" else ""} \
         ~{if (defined(parents) && select_first([parents])) then "-parents" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "breons/allsorts:0.1.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out_predictions = glob("predictions.csv")[0]
       File out_probabilities = glob("probabilities.csv")[0]
       File out_distributions = glob("distributions.png")[0]
       File out_waterfalls = glob("waterfalls.png")[0]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Allsorts
   doc: |
     usage: ALLSorts [-h] -samples SAMPLES [-labels LABELS]
                     [-destination DESTINATION] [-test] [-train]
                     [-model_dir MODEL_DIR] [-njobs NJOBS] [-cv CV] [-verbose]
                     [-comparison] [-force] [-parents]
     ALLSorts CLI

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: breons/allsorts:0.1.0

   inputs:
   - id: samples
     label: samples
     doc: |-
       (-s)  Path to samples (rows) x genes (columns) csv file representing a raw counts matrix. Note: hg19 only supported currently, use other references at own risk.
     type: File
     inputBinding:
       prefix: -samples
       separate: true
   - id: labels
     label: labels
     doc: |-
       (-l)  (Optional) Path to samples true labels. CSV with samples (rows) x [sample id, label] (cols). This will enable re-labelling mode. Note: labels must reflect naming conventions used within this tool. View the ALLSorts GitHub Wiki for further details.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -labels
       separate: true
   - id: destination
     label: destination
     doc: (-d)  Path to where you want the final report to be saved.
     type: string
     default: .
     inputBinding:
       prefix: -destination
       separate: true
   - id: verbose
     label: verbose
     doc: (-v) (flag, default=False) Verbose. Print stage progress.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -verbose
       separate: true
   - id: comparison
     label: comparison
     doc: Rebuild comparisons for labelled visualisations.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -comparison
       separate: true
   - id: force
     label: force
     doc: (-f) (flag, default=False) Force. Bypass warnings without user confirmation.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -force
       separate: true
   - id: parents
     label: parents
     doc: |-
       (-p) Include parent meta-subtypes in predictions. Note: This may remove previously unclassified samples.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -parents
       separate: true

   outputs:
   - id: out_predictions
     label: out_predictions
     type: File
     outputBinding:
       glob: predictions.csv
       loadContents: false
   - id: out_probabilities
     label: out_probabilities
     type: File
     outputBinding:
       glob: probabilities.csv
       loadContents: false
   - id: out_distributions
     label: out_distributions
     type: File
     outputBinding:
       glob: distributions.png
       loadContents: false
   - id: out_waterfalls
     label: out_waterfalls
     type: File
     outputBinding:
       glob: waterfalls.png
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - ALLSorts
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Allsorts


