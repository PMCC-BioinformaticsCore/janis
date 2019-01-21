
gatk4baserecalibrator
=====================
*bioinformatics*

Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php/>`_

Docstring
*********
First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
    The default covariates are read group, reported quality score, machine cycle, and nucleotide context.
    
    This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
    that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
    We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
    Since there is a large amount of data one can then calculate an empirical probability of error given the 
    particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
    table (of the several covariate values, num observations, num mismatches, empirical quality score).

*This page was automatically generated*
