
bgzip
=====
*bioinformatics*

Documentation
-------------

URL
******
`http://www.htslib.org/doc/bgzip.html <http://www.htslib.org/doc/bgzip.html/>`_

Docstring
*********
bgzip – Block compression/decompression utility

    Bgzip compresses files in a similar manner to, and compatible with, gzip(1). The file is compressed 
    into a series of small (less than 64K) 'BGZF' blocks. This allows indexes to be built against the 
    compressed file and used to retrieve portions of the data without having to decompress the entire file.

    If no files are specified on the command line, bgzip will compress (or decompress if the -d option is used) 
    standard input to standard output. If a file is specified, it will be compressed (or decompressed with -d). 
    If the -c option is used, the result will be written to standard output, otherwise when compressing bgzip 
    will write to a new file with a .gz suffix and remove the original. When decompressing the input file must 
    have a .gz suffix, which will be removed to make the output name. 
    Again after decompression completes the input file will be removed.

*This page was automatically generated*
