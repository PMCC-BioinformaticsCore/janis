
.. include:: bcftoolssort_v1.9

BCFTools: Sort
=============================

Description
-------------

Tool identifier: ``bcftoolssort``

Tool path: ``janis_bioinformatics.tools.bcftools.sort.sort_1_9 import BcfToolsSort_1_9``

Version: v1.9

Docker: ``michaelfranklin/bcftools:1.9``



Documentation
-------------

URL
******
*No URL to the documentation was provided*

Tool documentation
******************
About:   Sort VCF/BCF file.
Usage:   bcftools sort [OPTIONS] <FILE.vcf>

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ====================
name    type    prefix      position  documentation
======  ======  ========  ==========  ====================
vcf     VCF                        1  The VCF file to sort
======  ======  ========  ==========  ====================

Optional inputs
***************

==============  ==================  =============  ==========  =======================================================================================
name            type                prefix         position    documentation
==============  ==================  =============  ==========  =======================================================================================
outputFilename  Optional<Filename>  --output-file              (-o) output file name [stdout]
outputType      Optional<String>    --output-type              (-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
tempDir         Optional<String>    --temp-dir                 (-T) temporary files [/tmp/bcftools-sort.XXXXXX/]
==============  ==================  =============  ==========  =======================================================================================


Metadata
********

Author: **Unknown**


*BCFTools: Sort was last updated on 2019-07-11*.
*This page was automatically generated on 2019-07-24*.
