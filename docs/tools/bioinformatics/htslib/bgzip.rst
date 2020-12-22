:orphan:

BGZip
=============

``bgzip`` · *1 contributor · 2 versions*

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

       from janis_bioinformatics.tools.htslib.bgzip.bgzip_1_9 import BGZip_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bgzip_step",
           BGZip_1_9(
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

       file: file




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
:Container: quay.io/biocontainers/htslib:1.9--ha228f0b_7
:Authors: Michael Franklin
:Citations: None
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  =============  ===============
name    type           documentation
======  =============  ===============
out     Gzipped<File>
======  =============  ===============


Additional configuration (inputs)
---------------------------------

==============  ==================  ============  ==========  ========================================================================================================================================================================================================================================================
name            type                prefix          position  documentation
==============  ==================  ============  ==========  ========================================================================================================================================================================================================================================================
file            File                                     100  File to bgzip compress
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

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bgzip {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File file
       String? outputFilename
       Int? offset
       Boolean? stdout
       Boolean? decompress
       Boolean? force
       Boolean? help
       Boolean? index
       File? indexName
       Int? compress
       Boolean? reindex
       Boolean? rebgzip
       Int? size
       Int? threads
     }
     command <<<
       set -e
       bgzip \
         ~{if defined(offset) then ("--offset " + offset) else ''} \
         ~{if select_first([stdout, true]) then "--stdout" else ""} \
         ~{if (defined(decompress) && select_first([decompress])) then "--decompress" else ""} \
         ~{if (defined(force) && select_first([force])) then "--force" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if (defined(index) && select_first([index])) then "--index" else ""} \
         ~{if defined(indexName) then ("--index-name '" + indexName + "'") else ""} \
         ~{if defined(compress) then ("--compress " + compress) else ''} \
         ~{if (defined(reindex) && select_first([reindex])) then "--reindex" else ""} \
         ~{if (defined(rebgzip) && select_first([rebgzip])) then "--rebgzip" else ""} \
         ~{if defined(size) then ("--size " + size) else ''} \
         ~{if defined(threads) then ("--threads " + threads) else ''} \
         '~{file}' \
         > \
         '~{select_first([outputFilename, "~{basename(file)}.gz"])}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/htslib:1.9--ha228f0b_7"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(file)}.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: BGZip
   doc: |-
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/htslib:1.9--ha228f0b_7

   inputs:
   - id: file
     label: file
     doc: File to bgzip compress
     type: File
     inputBinding:
       position: 100
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.gz
     inputBinding:
       position: 102
       valueFrom: $(inputs.file.basename.basename).gz
   - id: offset
     label: offset
     doc: |-
       b: Decompress to standard output from virtual file position (0-based uncompressed offset). Implies -c and -d.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --offset
   - id: stdout
     label: stdout
     doc: 'c: Write to standard output, keep original files unchanged.'
     type: boolean
     default: true
     inputBinding:
       prefix: --stdout
   - id: decompress
     label: decompress
     doc: 'd: Decompress.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --decompress
   - id: force
     label: force
     doc: 'f: Overwrite files without asking.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force
   - id: help
     label: help
     doc: 'h: Displays a help message.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
   - id: index
     label: index
     doc: |-
       i: Create a BGZF index while compressing. Unless the -I option is used, this will have the name of the compressed file with .gzi appended to it.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --index
   - id: indexName
     label: indexName
     doc: '-I: Index file name.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --index-name
   - id: compress
     label: compress
     doc: |-
       l: Compression level to use when compressing. From 0 to 9, or -1 for the default level set by the compression library. [-1]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --compress
   - id: reindex
     label: reindex
     doc: 'r: Rebuild the index on an existing compressed file.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --reindex
   - id: rebgzip
     label: rebgzip
     doc: |-
       g: Try to use an existing index to create a compressed file with matching block offsets. Note that this assumes that the same compression library and level are in use as when making the original file. Don't use it unless you know what you're doing.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --rebgzip
   - id: size
     label: size
     doc: 's: Decompress INT bytes (uncompressed size) to standard output. Implies -c.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --size
   - id: threads
     label: threads
     doc: '@: Number of threads to use [1].'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: $(inputs.file.basename.basename).gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: bgzip
   arguments:
   - position: 101
     valueFrom: '>'
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bgzip


