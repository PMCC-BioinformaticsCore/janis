
.. include:: gatk4baserecalibrator_4.0.12.0

GATK4: Base Recalibrator
================================================

Description
-------------

Tool identifier: ``Gatk4BaseRecalibrator``

Tool path: ``janis_bioinformatics.tools.gatk4.baserecalibrator.latest import Gatk4BaseRecalibratorLatest``

Version: 4.0.12.0

Docker: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php>`_

Description
*********
First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
Since there is a large amount of data one can then calculate an empirical probability of error given the 
particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
table (of the several covariate values, num observations, num mismatches, empirical quality score).

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

==========  =================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name        type               prefix           position  documentation
==========  =================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bam         BamPair            -I                      6  BAM/SAM/CRAM file containing reads
knownSites  Array<vcf-gz-tbi>  --known-sites          28  **One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.** This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference, so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
reference   FastaWithDict      -R                      5  Reference sequence file
==========  =================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Optional inputs
***************

==============  ==================  ===========  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name            type                prefix         position  documentation
==============  ==================  ===========  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
tmpDir          Optional<String>    --tmp-dir                Temp directory to use.
outputFilename  Optional<Filename>  -O                    8  **The output recalibration table filename to create.** After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate. Use '/dev/stdout' to print to standard out.
intervals       Optional<bed>       --intervals              -L (BASE) One or more genomic intervals over which to operate
==============  ==================  ===========  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: Base Recalibrator was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-24*.
