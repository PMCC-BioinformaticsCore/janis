
.. include:: samtoolssort_1.9.0
.. include:: samtoolssort_1.7.0

SamTools: Sort
=============================

Description
-------------

Tool identifier: ``SamToolsSort``

Tool path: ``janis_bioinformatics.tools.samtools.sort.sort import SamToolsSort_1_9``

Version: 1.9.0

Container: ``quay.io/biocontainers/samtools:1.9--h8571acd_11``

Versions
*********

- 1.9.0 (current)
- `1.7.0 <samtoolssort_1.7.0.html>`_

Documentation
-------------

URL
******
`http://www.htslib.org/doc/samtools.html#DESCRIPTION <http://www.htslib.org/doc/samtools.html#DESCRIPTION>`_

Tool documentation
******************
Ensure SAMTOOLS.SORT is inheriting from parent metadata
    
---------------------------------------------------------------------------------------------------

Sort alignments by leftmost coordinates, or by read name when -n is used. An appropriate 
@HD-SO sort order header tag will be added or an existing one updated if necessary.

The sorted output is written to standard output by default, or to the specified file (out.bam) 
when -o is used. This command will also create temporary files tmpprefix.%d.bam as needed when 
the entire alignment data cannot fit into memory (as controlled via the -m option).

---------------------------------------------------------------------------------------------------

The following rules are used for ordering records.

If option -t is in use, records are first sorted by the value of the given alignment tag, and then 
by position or name (if using -n). For example, “-t RG” will make read group the primary sort key. 
The rules for ordering by tag are:

- Records that do not have the tag are sorted before ones that do.
- If the types of the tags are different, they will be sorted so that single character tags (type A) 
    come before array tags (type B), then string tags (types H and Z), then numeric tags (types f and i).
- Numeric tags (types f and i) are compared by value. Note that comparisons of floating-point values 
    are subject to issues of rounding and precision.
- String tags (types H and Z) are compared based on the binary contents of the tag using the C strcmp(3) function.
- Character tags (type A) are compared by binary character value.
- No attempt is made to compare tags of other types — notably type B array values will not be compared.

When the -n option is present, records are sorted by name. Names are compared so as to give a 
“natural” ordering — i.e. sections consisting of digits are compared numerically while all other 
sections are compared based on their binary representation. This means “a1” will come before 
“b1” and “a9” will come before “a10”. Records with the same name will be ordered according to 
the values of the READ1 and READ2 flags (see flags).

When the -n option is not present, reads are sorted by reference (according to the order of the 
@SQ header records), then by position in the reference, and then by the REVERSE flag.

*Note*

    Historically samtools sort also accepted a less flexible way of specifying the 
    final and temporary output filenames:
    
    |   samtools sort [-f] [-o] in.bam out.prefix
    
    This has now been removed. The previous out.prefix argument (and -f option, if any) 
    should be changed to an appropriate combination of -T PREFIX and -o FILE. The previous -o 
    option should be removed, as output defaults to standard output.

Outputs
-------
================  ====================  ==============================================================================================================================================================================================
name              type                  documentation
================  ====================  ==============================================================================================================================================================================================
out               BAM
temporaryOutputs  Optional<Array<BAM>>  By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or if output is to standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.
================  ====================  ==============================================================================================================================================================================================

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
bam     BAM                       10
======  ======  ========  ==========  ===============

Optional inputs
***************

====================  ==================  ========  ==========  ===========================================================================================================================================================================================================================================
name                  type                prefix      position  documentation
====================  ==================  ========  ==========  ===========================================================================================================================================================================================================================================
compression           Optional<Integer>   -l                    Set the desired compression level for the final output file, ranging from 0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best compression but slowest to write), similarly to gzip(1)'s compression level setting.
                                                                If -l is not used, the default compression level will apply.
maximumMemory         Optional<String>    -m                    Approximately the maximum required memory per thread, specified  either in bytes or with a K, M, or G suffix [768 MiB]. To prevent sort from creating a huge number of temporary files, it enforces a minimum value of 1M for this setting.
sortByReadNames       Optional<Boolean>   -n                    Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
outputType            Optional<String>    -O                    Write the final output as sam, bam, or cram. By default, samtools tries to select a format based on the -o filename extension; if output is to standard output or no format can be deduced, bam is selected.
temporaryFilesPrefix  Optional<String>    -T                    Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique to this invocation of the sort command.
                                                                By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or if output is to standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.
threads               Optional<Integer>   -@                    Set number of sorting and compression threads. By default, operation is single-threaded.
outputFilename        Optional<Filename>  -o                 5  Output to FILE [stdout].
====================  ==================  ========  ==========  ===========================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*SamTools: Sort was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-30*.
