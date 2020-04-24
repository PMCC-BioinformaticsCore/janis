:orphan:

GATK4: Print Reads
====================================

*1 contributor · 4 versions*


Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
A common use case is to subset reads by genomic interval using the -L argument. 
Note when applying genomic intervals, the tool is literal and does not retain mates 
of paired-end reads outside of the interval, if any. Data with missing mates will fail 
ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. 
If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

By default, PrintReads applies the WellformedReadFilter at the engine level. 
What this means is that the tool does not print reads that fail the WellformedReadFilter filter. 
You can similarly apply other engine-level filters to remove specific types of reads 
with the --read-filter argument. See documentation category 'Read Filters' for a list of
 available filters. To keep reads that do not pass the WellformedReadFilter, either 
 disable the filter with --disable-read-filter or disable all default filters with 
 ``--disable-tool-default-read-filters``.

The reference is strictly required when handling CRAM files.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.printreads.versions import Gatk4PrintReads_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4printreads_step",
           Gatk4PrintReads_4_1_3(
               bam=None,
           )
       )
       wf.output("out", source=gatk4printreads_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4PrintReads:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4PrintReads > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run Gatk4PrintReads with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4PrintReads





Information
------------


:ID: ``Gatk4PrintReads``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ===============
name            type                prefix    position    documentation
==============  ==================  ========  ==========  ===============
bam             BAM
outputFilename  Optional<Filename>
==============  ==================  ========  ==========  ===============
