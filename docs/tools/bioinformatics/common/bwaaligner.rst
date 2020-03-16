:orphan:

Align and sort reads
=================================

*0 contributors · 1 version*

Align sorted bam with this subworkflow consisting of BWA Mem + SamTools + Gatk4SortSam

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.bwaaligner import BwaAligner

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bwaaligner_step",
           BwaAligner(
               sample_name=None,
               reference=None,
               fastq=None,
           )
       )
       wf.output("out", source=bwaaligner_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for BwaAligner:

.. code-block:: bash

   # user inputs
   janis inputs BwaAligner > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastq:
       - fastq_0.fastq.gz
       - fastq_1.fastq.gz
       reference: reference.fasta
       sample_name: <value>




5. Run BwaAligner with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       BwaAligner





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``BwaAligner``
:URL: *No URL to the documentation was provided*
:Versions: 1.0.0
:Authors: 
:Citations: 
:Created: 2018-12-24
:Updated: None



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Embedded Tools
***************

=======================  =================================
Cutadapt                 ``cutadapt/2.6``
Bwa mem + Samtools View  ``BwaMemSamtoolsView/0.7.17|1.9``
GATK4: SortSAM           ``Gatk4SortSam/4.1.3.0``
=======================  =================================



Additional configuration (inputs)
---------------------------------

=============================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================
name                           type                     documentation
=============================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================
sample_name                    String
reference                      FastaWithIndexes
fastq                          FastqGzPair
cutadapt_adapter               Optional<Array<String>>
cutadapt_removeMiddle3Adapter  Optional<Array<String>>
cutadapt_front                 Optional<String>
cutadapt_removeMiddle5Adapter  Optional<String>
cutadapt_qualityCutoff         Optional<Integer>        (]3'CUTOFF, ]3'CUTOFF, -q)  Trim low-quality bases from 5' and/or 3' ends of each read before adapter removal. Applied to both reads if data is paired. If one value is given, only the 3' end is trimmed. If two comma-separated cutoffs are given, the 5' end is trimmed with the first cutoff, the 3' end with the second.
cutadapt_minimumLength         Optional<Integer>        (-m)  Discard reads shorter than LEN. Default: 0
bwamem_markShorterSplits       Optional<Boolean>        Mark shorter split hits as secondary (for Picard compatibility).
sortsam_sortOrder              Optional<String>         The --SORT_ORDER argument is an enumerated type (SortOrder), which can have one of the following values: [unsorted, queryname, coordinate, duplicate, unknown]
sortsam_createIndex            Optional<Boolean>        Whether to create a BAM index when writing a coordinate-sorted BAM file.
sortsam_validationStringency   Optional<String>         Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
sortsam_maxRecordsInRam        Optional<Integer>        When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
sortsam_tmpDir                 Optional<String>         Undocumented option
=============================  =======================  ================================================================================================================================================================================================================================================================================================================================================================================================


