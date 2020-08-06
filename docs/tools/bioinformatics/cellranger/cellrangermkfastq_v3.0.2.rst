:orphan:

CellRanger mkfastq
======================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.cellranger.mkfastq.versions import CellRangerMkfastq_3_0_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "cellrangermkfastq_step",
           CellRangerMkfastq_3_0_2(
               run=None,
           )
       )
       wf.output("out", source=cellrangermkfastq_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for CellRangerMkfastq:

.. code-block:: bash

   # user inputs
   janis inputs CellRangerMkfastq > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       run: null




5. Run CellRangerMkfastq with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       CellRangerMkfastq





Information
------------

:ID: ``CellRangerMkfastq``
:URL: *No URL to the documentation was provided*
:Versions: v3.0.2
:Container: fbrundu/cellranger:v3.0.2
:Authors: 
:Citations: None
:Created: None
:Updated: None


Outputs
-----------

======  =========  ===============
name    type       documentation
======  =========  ===============
out     Directory
======  =========  ===============


Additional configuration (inputs)
---------------------------------

==================  =======================  =====================  ==========  ================================================================================================================================================================================================================================================================================
name                type                     prefix                 position    documentation
==================  =======================  =====================  ==========  ================================================================================================================================================================================================================================================================================
run                 Directory                --run=                             Path of Illumina BCL run folder.
id                  Optional<String>         --id=                              Name of the folder created by mkfastq. If not supplied, will default to the name of the flowcell referred to by the --run argument.
outputFoldername    Optional<Filename>       --output-dir=                      Same as in bcl2fastq. Folder where FASTQs, reports and stats will be generated.
csv                 Optional<csv>            --csv=                             Apparently the same as `sampleSheet`. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
sampleSheet         Optional<File>           --sample-sheet=                    (--samplesheet= | --csv=) Path to the sample sheet. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
ignoreDualIndex     Optional<Boolean>        --ignore-dual-index                On a dual-indexed flowcell, ignore the second sample index, if the second sample index was not used for the 10x sample.
qc                  Optional<Boolean>        --qc                               Calculate both sequencing and 10x-specific metrics, including per-sample barcode matching rate. Will not be performed unless this flag is specified.
lanes               Optional<Array<String>>  --lanes=                           Comma-delimited series of lanes to demultiplex. Shortcut for the --tiles argument.
useBasesMask        Optional<String>         --use-bases-mask=                  Same as bcl2fastq; override the read lengths as specified in RunInfo.xml. See Illumina bcl2fastq documentation for more information.
deleteUndetermined  Optional<Boolean>        --delete-undetermined              Delete the Undetermined FASTQ files left by bcl2fastq.  Useful if your sample sheet is only expected to match a subset of the flowcell.
project             Optional<String>         --project=                         Custom project name, to override the samplesheet or to use in conjunction with the --csv argument.
localcores          Optional<Integer>        --localcores=                      Set max cores the pipeline may request at one time. Only applies when --jobmode=local.
localmem            Optional<Float>          --localmem=                        Set max GB the pipeline may request at one time. Only applies when --jobmode=local.
nopreflight         Optional<Boolean>        --nopreflight                      Skip preflight checks.
==================  =======================  =====================  ==========  ================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task CellRangerMkfastq {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Directory run
       String? id
       String? outputFoldername
       File? csv
       File? sampleSheet
       Boolean? ignoreDualIndex
       Boolean? qc
       Array[String]? lanes
       String? useBasesMask
       Boolean? deleteUndetermined
       String? project
       Int? localcores
       Float? localmem
       Boolean? nopreflight
     }
     command <<<
       set -e
       cellranger mkfastq \
         --run='~{run}' \
         ~{if defined(id) then ("--id='" + id + "'") else ""} \
         --output-dir='~{select_first([outputFoldername, "generated"])}' \
         ~{if defined(csv) then ("--csv='" + csv + "'") else ""} \
         ~{if defined(sampleSheet) then ("--sample-sheet='" + sampleSheet + "'") else ""} \
         ~{if defined(ignoreDualIndex) then "--ignore-dual-index" else ""} \
         ~{if defined(qc) then "--qc" else ""} \
         ~{if (defined(lanes) && length(select_first([lanes])) > 0) then "--lanes='" + sep("','", select_first([lanes])) + "'" else ""} \
         ~{if defined(useBasesMask) then ("--use-bases-mask='" + useBasesMask + "'") else ""} \
         ~{if defined(deleteUndetermined) then "--delete-undetermined" else ""} \
         ~{if defined(project) then ("--project='" + project + "'") else ""} \
         ~{if defined(select_first([localcores, select_first([runtime_cpu, 1])])) then ("--localcores=" + select_first([localcores, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(localmem) then ("--localmem=" + localmem) else ''} \
         ~{if defined(nopreflight) then "--nopreflight" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "fbrundu/cellranger:v3.0.2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       Directory out = select_first([outputFoldername, "generated"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.0
   label: CellRanger mkfastq

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: fbrundu/cellranger:v3.0.2

   inputs:
   - id: run
     label: run
     doc: Path of Illumina BCL run folder.
     type: Directory
     inputBinding:
       prefix: --run=
       separate: false
   - id: id
     label: id
     doc: |-
       Name of the folder created by mkfastq. If not supplied, will default to the name of the flowcell referred to by the --run argument.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --id=
       separate: false
   - id: outputFoldername
     label: outputFoldername
     doc: Same as in bcl2fastq. Folder where FASTQs, reports and stats will be generated.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --output-dir=
       separate: false
   - id: csv
     label: csv
     doc: |-
       Apparently the same as `sampleSheet`. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --csv=
       separate: false
   - id: sampleSheet
     label: sampleSheet
     doc: |-
       (--samplesheet= | --csv=) Path to the sample sheet. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --sample-sheet=
       separate: false
   - id: ignoreDualIndex
     label: ignoreDualIndex
     doc: |-
       On a dual-indexed flowcell, ignore the second sample index, if the second sample index was not used for the 10x sample.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ignore-dual-index
       separate: true
   - id: qc
     label: qc
     doc: |-
       Calculate both sequencing and 10x-specific metrics, including per-sample barcode matching rate. Will not be performed unless this flag is specified.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --qc
       separate: true
   - id: lanes
     label: lanes
     doc: |-
       Comma-delimited series of lanes to demultiplex. Shortcut for the --tiles argument.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --lanes=
       separate: false
       itemSeparator: ','
   - id: useBasesMask
     label: useBasesMask
     doc: |-
       Same as bcl2fastq; override the read lengths as specified in RunInfo.xml. See Illumina bcl2fastq documentation for more information.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --use-bases-mask=
       separate: false
   - id: deleteUndetermined
     label: deleteUndetermined
     doc: |-
       Delete the Undetermined FASTQ files left by bcl2fastq.  Useful if your sample sheet is only expected to match a subset of the flowcell.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --delete-undetermined
       separate: true
   - id: project
     label: project
     doc: |-
       Custom project name, to override the samplesheet or to use in conjunction with the --csv argument.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --project=
       separate: false
   - id: localcores
     label: localcores
     doc: |-
       Set max cores the pipeline may request at one time. Only applies when --jobmode=local.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --localcores=
       valueFrom: $([inputs.runtime_cpu, 1].filter(function (inner) { return inner !=
         null })[0])
       separate: false
   - id: localmem
     label: localmem
     doc: |-
       Set max GB the pipeline may request at one time. Only applies when --jobmode=local.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --localmem=
       valueFrom: |-
         $([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0])
       separate: false
   - id: nopreflight
     label: nopreflight
     doc: Skip preflight checks.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --nopreflight
       separate: true

   outputs:
   - id: out
     label: out
     type: Directory
     outputBinding:
       glob: generated
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - cellranger
   - mkfastq
   arguments: []
   id: CellRangerMkfastq


