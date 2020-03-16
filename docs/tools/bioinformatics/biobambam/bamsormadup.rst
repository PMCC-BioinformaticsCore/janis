:orphan:

BamSorMaDup
=========================

*1 contributor · 1 version*

:ID: ``bamsormadup``
:Python: ``janis_bioinformatics.tools.biobambam.bamsormadup.versions import BamSorMaDup_2_0_87``
:Versions: 2.0.87
:Container: quay.io/biocontainers/biobambam:2.0.87--1
:Authors: Matthias De Smet (@mattdsm)
:Citations: None
:Created: 2020-02-26
:Updated: 2020-02-26
:Required inputs:
   - ``alignedReads: BAM``
:Outputs: 
   - ``out: stdout<BAM>``

   - ``metrics: File``

Documentation
-------------

URL: `https://gitlab.com/german.tischler/biobambam2 <https://gitlab.com/german.tischler/biobambam2>`_

bamsormadup: parallel sorting and duplicate marking

------

Arguments
----------

===========  =============  ==========  ==============================================
value        prefix           position  documentation
===========  =============  ==========  ==============================================
metrics.txt  M=                      0  file containing metrics from duplicate removal
bam          inputformat=            0  input data format
bam          outputFormat=           0  output data format
===========  =============  ==========  ==============================================

Additional configuration (inputs)
---------------------------------

==============  ==================  ===============  ==========  =========================================================================================================
name            type                prefix             position  documentation
==============  ==================  ===============  ==========  =========================================================================================================
alignedReads    BAM                                         200
outputFilename  Optional<Filename>
level           Optional<Integer>   level=                       compression settings for output bam file (-1=zlib default,0=uncompressed,1=fast,9=best)
tempLevel       Optional<Integer>   templevel=                   compression settings for temporary bam files (-1=zlib default,0=uncompressed,1=fast,9=best)
threads         Optional<Integer>   threads=                     Number of threads. (default = 1)
sortOrder       Optional<String>    SO=                          output sort order(coordinate by default)
optMinPixelDif  Optional<Integer>   optminpixeldif=              pixel difference threshold for optical duplicates (patterned flowcell: 12000, unpatterned flowcell: 2500)
==============  ==================  ===============  ==========  =========================================================================================================

