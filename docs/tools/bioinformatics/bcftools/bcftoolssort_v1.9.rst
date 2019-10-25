:orphan:

BCFTools: Sort
=============================

0 contributors · 1 version

:ID: ``bcftoolssort``
:Python: ``janis_bioinformatics.tools.bcftools.sort.versions import BcfToolsSort_1_9``
:Versions: v1.9
:Container: michaelfranklin/bcftools:1.9
:Authors: 
:Citations: None
:Created: 2019-05-09
:Updated: 2019-07-11
:Required inputs:
   - ``vcf: CompressedVCF``
:Outputs: 
   - ``out: CompressedVCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

About:   Sort VCF/BCF file.
Usage:   bcftools sort [OPTIONS] <FILE.vcf>

------

Additional configuration (inputs)
---------------------------------

==============  ==================  =======================================================================================
name            type                documentation
==============  ==================  =======================================================================================
vcf             CompressedVCF       The VCF file to sort
outputFilename  Optional<Filename>  (-o) output file name [stdout]
outputType      Optional<String>    (-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
tempDir         Optional<String>    (-T) temporary files [/tmp/bcftools-sort.XXXXXX/]
==============  ==================  =======================================================================================

