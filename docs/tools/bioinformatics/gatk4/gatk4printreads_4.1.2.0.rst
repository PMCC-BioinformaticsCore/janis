:orphan:


GATK4: Print Reads
====================================

Description
-------------

Tool identifier: ``Gatk4PrintReads``

Tool path: ``janis_bioinformatics.tools.gatk4.printreads.versions import Gatk4PrintReads_4_1_2``

Version: 4.1.2.0

Container: ``broadinstitute/gatk:4.1.2.0``

Versions
*********

- `4.1.3.0 <gatk4printreads_4.1.3.0.html>`_
- 4.1.2.0 (current)
- `4.0.12.0 <gatk4printreads_4.0.12.0.html>`_

Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php>`_

Tool documentation
******************

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

Outputs
-------
======  =======  ===============
name    type     documentation
======  =======  ===============
out     BamPair
======  =======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix    position    documentation
======  ======  ========  ==========  ===============
bam     BAM
======  ======  ========  ==========  ===============

Optional inputs
***************

==============  ==================  ========  ==========  ===============
name            type                prefix    position    documentation
==============  ==================  ========  ==========  ===============
outputFilename  Optional<Filename>
==============  ==================  ========  ==========  ===============


Metadata
********

Author: Michael Franklin


*GATK4: Print Reads was last updated on 2019-01-24*.
*This page was automatically generated on 2019-09-10*.