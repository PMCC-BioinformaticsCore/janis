:orphan:

CellRanger mkfastq
======================================

0 contributors · 1 version

:ID: ``CellRangerMkfastq``
:Python: ``janis_bioinformatics.tools.cellranger.mkfastq.versions import CellRangerMkfastq_3_0_2``
:Versions: v3.0.2
:Container: fbrundu/cellranger:v3.0.2
:Authors: 
:Citations: None
:Created: None
:Updated: None
:Required inputs:
   - ``run: Directory``
:Outputs: 
   - ``out: Directory``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

------

Additional configuration (inputs)
---------------------------------

==================  =======================  ================================================================================================================================================================================================================================================================================
name                type                     documentation
==================  =======================  ================================================================================================================================================================================================================================================================================
run                 Directory                Path of Illumina BCL run folder.
id                  Optional<String>         Name of the folder created by mkfastq. If not supplied, will default to the name of the flowcell referred to by the --run argument.
outputFoldername    Optional<Filename>       Same as in bcl2fastq. Folder where FASTQs, reports and stats will be generated.
csv                 Optional<csv>            Apparently the same as `sampleSheet`. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
sampleSheet         Optional<File>           (--samplesheet= | --csv=) Path to the sample sheet. The sample sheet can either be a simple CSV with lane, sample and index columns, or an Illumina Experiment Manager-compatible sample sheet.  Sample sheet indexes can refer to 10x sample index set names (e.g., SI-GA-A12).
ignoreDualIndex     Optional<Boolean>        On a dual-indexed flowcell, ignore the second sample index, if the second sample index was not used for the 10x sample.
qc                  Optional<Boolean>        Calculate both sequencing and 10x-specific metrics, including per-sample barcode matching rate. Will not be performed unless this flag is specified.
lanes               Optional<Array<String>>  Comma-delimited series of lanes to demultiplex. Shortcut for the --tiles argument.
useBasesMask        Optional<String>         Same as bcl2fastq; override the read lengths as specified in RunInfo.xml. See Illumina bcl2fastq documentation for more information.
deleteUndetermined  Optional<Boolean>        Delete the Undetermined FASTQ files left by bcl2fastq.  Useful if your sample sheet is only expected to match a subset of the flowcell.
project             Optional<String>         Custom project name, to override the samplesheet or to use in conjunction with the --csv argument.
localcores          Optional<Integer>        Set max cores the pipeline may request at one time. Only applies when --jobmode=local.
localmem            Optional<Float>          Set max GB the pipeline may request at one time. Only applies when --jobmode=local.
nopreflight         Optional<Boolean>        Skip preflight checks.
==================  =======================  ================================================================================================================================================================================================================================================================================

