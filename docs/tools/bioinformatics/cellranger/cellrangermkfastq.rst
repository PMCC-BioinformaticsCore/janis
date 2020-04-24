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
