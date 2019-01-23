
gatkprintreads
==============
*bioinformatics*

Documentation
-------------

URL
******
*No URL to the documentation was provided*: `contribute one <https://github.com/illusional>`_

Docstring
*********
Write reads from SAM format file (SAM/BAM/CRAM) that pass criteria to a new file.
A common use case is to subset reads by genomic interval using the -L argument. 
Note when applying genomic intervals, the tool is literal and does not retain mates 
of paired-end reads outside of the interval, if any. Data with missing mates will fail 
ValidateSamFile validation with MATE_NOT_FOUND, but certain tools may still analyze the data. 
If needed, to rescue such mates, use either FilterSamReads or ExtractOriginalAlignmentRecordsByNameSpark.

By default, PrintReads applies the WellformedReadFilter at the engine level. What this means is that 
the tool does not print reads that fail the WellformedReadFilter filter. You can similarly apply 
other engine-level filters to remove specific types of reads with the --read-filter argument. 
See documentation category 'Read Filters' for a list of available filters. 
To keep reads that do not pass the WellformedReadFilter, either disable the filter 
with --disable-read-filter or disable all default filters with --disable-tool-default-read-filters.

The reference is strictly required when handling CRAM files.

Documentation: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_PrintReads.php

Outputs
-------
============  =======  =========================
name          type     documentation
============  =======  =========================
pairedOutput  BamPair  Write output to this file
============  =======  =========================

Inputs
------
======================  =======================  ===================  ==========  ============================================================================
name                    type                     prefix                 position  documentation
======================  =======================  ===================  ==========  ============================================================================
showFullBamList         Optional<Boolean>        --showFullBamList                Emit list of input BAM/CRAM files to log
inputBam                BamPair                  -I                            6  BAM/SAM/CRAM file containing reads
reference               FastaWithDict            -R                            5
input_baseRecalibrator  File                     -BQSR                         7  the recalibration table produced by BaseRecalibration
bedFile                 bed                      -L                           15
sample_file             Optional<Array<File>>                                 11
platform                Optional<String>         --platform                   13  Exclude all reads with this platform from the output
number                  Optional<String>         --number                     13  Exclude all reads with this platform from the output
simplify                Optional<Boolean>        --simplify                    9  Erase all extra attributes in the read but keep the read group information
readGroup               Optional<String>         --readGroup                  12  Exclude all reads with this read group from the output
sample_name             Optional<Array<String>>                               10  Sample name to be included in the analysis. Can be specified multiple times.
outputFilename          Optional<String>         -o                            8  name of the output file from indelRealigner
java_arg                Optional<String>                                       1
threads                 Optional<Integer>        -nct                         14
downsamplingType        Optional<String>         --downsampling_type          16
======================  =======================  ===================  ==========  ============================================================================


*This page was automatically generated*
