:orphan:


BGZip
=============

Description
-------------

Tool identifier: ``bgzip``

Tool path: ``janis_bioinformatics.tools.htslib.bgzip.bgzip_1_9 import BGZip_1_9``

Version: 1.9

Docker: ``quay.io/biocontainers/htslib:1.9--ha228f0b_7``

Versions
*********

- 1.9 (current)
- `1.2.1 <bgzip_1.2.1.html>`_

Documentation
-------------

URL
******
`http://www.htslib.org/doc/bgzip.html <http://www.htslib.org/doc/bgzip.html>`_

Tool documentation
******************
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

Outputs
-------
======  ======  ===============
name    type    documentation
======  ======  ===============
out     Stdout
======  ======  ===============

Inputs
------
Find the inputs below

Required inputs
***************

======  ======  ========  ==========  ======================
name    type    prefix      position  documentation
======  ======  ========  ==========  ======================
file    VCF                      100  File to bgzip compress
======  ======  ========  ==========  ======================

Optional inputs
***************

==========  =================  ============  ==========  ========================================================================================================================================================================================================================================================
name        type               prefix        position    documentation
==========  =================  ============  ==========  ========================================================================================================================================================================================================================================================
offset      Optional<Integer>  --offset                  b: Decompress to standard output from virtual file position (0-based uncompressed offset). Implies -c and -d.
stdout      Optional<Boolean>  --stdout                  c: Write to standard output, keep original files unchanged.
decompress  Optional<Boolean>  --decompress              d: Decompress.
force       Optional<Boolean>  --force                   f: Overwrite files without asking.
help        Optional<Boolean>  --help                    h: Displays a help message.
index       Optional<Boolean>  --index                   i: Create a BGZF index while compressing. Unless the -I option is used, this will have the name of the compressed file with .gzi appended to it.
indexName   Optional<File>     --index-name              -I: Index file name.
compress    Optional<Integer>  --compress                l: Compression level to use when compressing. From 0 to 9, or -1 for the default level set by the compression library. [-1]
reindex     Optional<Boolean>  --reindex                 r: Rebuild the index on an existing compressed file.
rebgzip     Optional<Boolean>  --rebgzip                 g: Try to use an existing index to create a compressed file with matching block offsets. Note that this assumes that the same compression library and level are in use as when making the original file. Don't use it unless you know what you're doing.
size        Optional<Integer>  --size                    s: Decompress INT bytes (uncompressed size) to standard output. Implies -c.
threads     Optional<Integer>  --threads                 @: Number of threads to use [1].
==========  =================  ============  ==========  ========================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*BGZip was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-29*.
