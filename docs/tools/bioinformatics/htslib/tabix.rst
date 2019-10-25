:orphan:

Tabix
=============

1 contributor · 2 versions

:ID: ``tabix``
:Python: ``janis_bioinformatics.tools.htslib.tabix.tabix_1_9 import Tabix_1_9``
:Versions: 1.9, 1.2.1
:Container: quay.io/biocontainers/htslib:1.9--ha228f0b_7
:Authors: Michael Franklin
:Citations: None
:Created: 2018-12-24
:Updated: 2019-01-24
:Required inputs:
   - ``file: CompressedVCF``
:Outputs: 
   - ``out: CompressedIndexedVCF``

Documentation
-------------

URL: `http://www.htslib.org/doc/tabix.html <http://www.htslib.org/doc/tabix.html>`_

tabix – Generic indexer for TAB-delimited genome position files

Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
and compressed by bgzip which has a gzip(1) like interface.

After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
"chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

Fast data retrieval also works over network if URI is given as a file name and in this case the 
index file will be downloaded if it is not present locally.

------

Additional configuration (inputs)
---------------------------------

===========  =================  ==============================================================================================================================================================================================================================================================================================================
name         type               documentation
===========  =================  ==============================================================================================================================================================================================================================================================================================================
file         CompressedVCF      File from which to create the index. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
preset       Optional<String>   -p: Input format for indexing. Valid values are: gff, bed, sam, vcf. This option should not be applied together with any of -s, -b, -e, -c and -0; it is not used for data retrieval because this setting is stored in the index file. [gff]
zeroBased    Optional<Boolean>  -0: Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based.
begin        Optional<Integer>  -b: Column of start chromosomal position. [4]
comment      Optional<String>   -c: Skip lines started with character CHAR. [#]
csi          Optional<Boolean>  -C: Produce CSI format index instead of classical tabix or BAI style indices.
end          Optional<Integer>  -e: Column of end chromosomal position. The end column can be the same as the start column. [5]
force        Optional<Boolean>  -f: Force to overwrite the index file if it is present.
minShift     Optional<Integer>  -m: set minimal interval size for CSI indices to 2^INT [14]
sequence     Optional<Integer>  -s: Column of sequence name. Option -s, -b, -e, -S, -c and -0 are all stored in the index file and thus not used in data retrieval. [1]
skipLines    Optional<Integer>  -S: Skip first INT lines in the data file. [0]
printHeader  Optional<Boolean>  -h: Print also the header/meta lines.
onlyHeader   Optional<Boolean>  -H: Print only the header/meta lines.
listChroms   Optional<Boolean>  -l: List the sequence names stored in the index file.
reheader     Optional<File>     -r: Replace the header with the content of FILE
regions      Optional<File>     -R: Restrict to regions listed in the FILE. The FILE can be BED file (requires .bed, .bed.gz, .bed.bgz file name extension) or a TAB-delimited file with CHROM, POS, and, optionally, POS_TO columns, where positions are 1-based and inclusive. When this option is in use, the input file may not be sorted.
targets      Optional<File>     -T: Similar to -R but the entire input will be read sequentially and regions not listed in FILE will be skipped
===========  =================  ==============================================================================================================================================================================================================================================================================================================

