:orphan:


Tabix
=============

Description
-------------

Tool identifier: ``tabix``

Tool path: ``janis_bioinformatics.tools.htslib.tabix.latest import TabixLatest``

Version: 1.2.1

Docker: ``biodckrdev/htslib:1.2.1``

Versions
*********

- `1.9 <tabix_1.9.html>`_
- 1.2.1 (current)

Documentation
-------------

URL
******
`http://www.htslib.org/doc/tabix.html <http://www.htslib.org/doc/tabix.html>`_

Tool documentation
******************
tabix – Generic indexer for TAB-delimited genome position files

Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
and compressed by bgzip which has a gzip(1) like interface.

After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
"chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

Fast data retrieval also works over network if URI is given as a file name and in this case the 
index file will be downloaded if it is not present locally.

Outputs
-------
======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     vcf-gz-tbi
======  ==========  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  =================  ========  ==========  ============================================================================================================================================
name    type               prefix      position  documentation
======  =================  ========  ==========  ============================================================================================================================================
file    compressed-vcf-gz                     8  File from which to create the index. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
======  =================  ========  ==========  ============================================================================================================================================

Optional inputs
***************

===========  =================  ==============  ==========  ==============================================================================================================================================================================================================================================================================================================
name         type               prefix            position  documentation
===========  =================  ==============  ==========  ==============================================================================================================================================================================================================================================================================================================
preset       Optional<String>   --preset                 2  -p: Input format for indexing. Valid values are: gff, bed, sam, vcf. This option should not be applied together with any of -s, -b, -e, -c and -0; it is not used for data retrieval because this setting is stored in the index file. [gff]
zeroBased    Optional<Boolean>  --zero-based             1  -0: Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based.
begin        Optional<Integer>  --begin                  4  -b: Column of start chromosomal position. [4]
comment      Optional<String>   --comment                7  -c: Skip lines started with character CHAR. [#]
csi          Optional<Boolean>  --csi                    1  -C: Produce CSI format index instead of classical tabix or BAI style indices.
end          Optional<Integer>  --end                    5  -e: Column of end chromosomal position. The end column can be the same as the start column. [5]
force        Optional<Boolean>  --force                  1  -f: Force to overwrite the index file if it is present.
minShift     Optional<Integer>  --min-shift              1  -m: set minimal interval size for CSI indices to 2^INT [14]
sequence     Optional<Integer>  --sequence               3  -s: Column of sequence name. Option -s, -b, -e, -S, -c and -0 are all stored in the index file and thus not used in data retrieval. [1]
skipLines    Optional<Integer>  --skip-lines             6  -S: Skip first INT lines in the data file. [0]
printHeader  Optional<Boolean>  --print-header           1  -h: Print also the header/meta lines.
onlyHeader   Optional<Boolean>  --only-header            1  -H: Print only the header/meta lines.
listChroms   Optional<Boolean>  --list-chroms            1  -l: List the sequence names stored in the index file.
reheader     Optional<File>     --reheader               1  -r: Replace the header with the content of FILE
regions      Optional<File>     --regions               11  -R: Restrict to regions listed in the FILE. The FILE can be BED file (requires .bed, .bed.gz, .bed.bgz file name extension) or a TAB-delimited file with CHROM, POS, and, optionally, POS_TO columns, where positions are 1-based and inclusive. When this option is in use, the input file may not be sorted.
targets      Optional<File>     --targets               11  -T: Similar to -R but the entire input will be read sequentially and regions not listed in FILE will be skipped
===========  =================  ==============  ==========  ==============================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*Tabix was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-26*.
