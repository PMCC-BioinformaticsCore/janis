:orphan:

BGZip
=============

*1 contributor · 2 versions*

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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.htslib.bgzip.bgzip_1_2_1 import BGZip_1_2_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bgzip_step",
           BGZip_1_2_1(
               file=None,
           )
       )
       wf.output("out", source=bgzip_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bgzip:

.. code-block:: bash

   # user inputs
   janis inputs bgzip > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       file: file.vcf




5. Run bgzip with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bgzip





Information
------------


:ID: ``bgzip``
:URL: `http://www.htslib.org/doc/bgzip.html <http://www.htslib.org/doc/bgzip.html>`_
:Versions: 1.9, 1.2.1
:Container: biodckrdev/htslib:1.2.1
:Authors: Michael Franklin
:Citations: None
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  =============  ===============
name    type           documentation
======  =============  ===============
out     CompressedVCF
======  =============  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ============  ==========  ========================================================================================================================================================================================================================================================
name            type                prefix          position  documentation
==============  ==================  ============  ==========  ========================================================================================================================================================================================================================================================
file            VCF                                      100  File to bgzip compress
outputFilename  Optional<Filename>                       102
offset          Optional<Integer>   --offset                  b: Decompress to standard output from virtual file position (0-based uncompressed offset). Implies -c and -d.
stdout          Optional<Boolean>   --stdout                  c: Write to standard output, keep original files unchanged.
decompress      Optional<Boolean>   --decompress              d: Decompress.
force           Optional<Boolean>   --force                   f: Overwrite files without asking.
help            Optional<Boolean>   --help                    h: Displays a help message.
index           Optional<Boolean>   --index                   i: Create a BGZF index while compressing. Unless the -I option is used, this will have the name of the compressed file with .gzi appended to it.
indexName       Optional<File>      --index-name              -I: Index file name.
compress        Optional<Integer>   --compress                l: Compression level to use when compressing. From 0 to 9, or -1 for the default level set by the compression library. [-1]
reindex         Optional<Boolean>   --reindex                 r: Rebuild the index on an existing compressed file.
rebgzip         Optional<Boolean>   --rebgzip                 g: Try to use an existing index to create a compressed file with matching block offsets. Note that this assumes that the same compression library and level are in use as when making the original file. Don't use it unless you know what you're doing.
size            Optional<Integer>   --size                    s: Decompress INT bytes (uncompressed size) to standard output. Implies -c.
threads         Optional<Integer>   --threads                 @: Number of threads to use [1].
==============  ==================  ============  ==========  ========================================================================================================================================================================================================================================================
