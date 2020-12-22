:orphan:

SamTools: Sort
=============================

``SamToolsSort`` · *1 contributor · 2 versions*

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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.sort.sort import SamToolsSort_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolssort_step",
           SamToolsSort_1_7(
               bam=None,
           )
       )
       wf.output("out", source=samtoolssort_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsSort:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsSort > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run SamToolsSort with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsSort





Information
------------

:ID: ``SamToolsSort``
:URL: `http://www.htslib.org/doc/samtools.html#DESCRIPTION <http://www.htslib.org/doc/samtools.html#DESCRIPTION>`_
:Versions: 1.9.0, 1.7.0
:Container: biocontainers/samtools:v1.7.0_cv3
:Authors: Michael Franklin
:Citations: None
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     BAM
======  ======  ===============


Additional configuration (inputs)
---------------------------------

====================  ==================  ========  ==========  ===========================================================================================================================================================================================================================================
name                  type                prefix      position  documentation
====================  ==================  ========  ==========  ===========================================================================================================================================================================================================================================
bam                   BAM                                   10
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

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task SamToolsSort {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Int? compression
       String? maximumMemory
       Boolean? sortByReadNames
       String? outputType
       String? temporaryFilesPrefix
       Int? threads
       File bam
       String? outputFilename
     }
     command <<<
       set -e
       samtools sort \
         ~{if defined(compression) then ("-l " + compression) else ''} \
         ~{if defined(maximumMemory) then ("-m '" + maximumMemory + "'") else ""} \
         ~{if (defined(sortByReadNames) && select_first([sortByReadNames])) then "-n" else ""} \
         ~{if defined(outputType) then ("-O '" + outputType + "'") else ""} \
         ~{if defined(temporaryFilesPrefix) then ("-T '" + temporaryFilesPrefix + "'") else ""} \
         ~{if defined(threads) then ("-@ " + threads) else ''} \
         -o '~{select_first([outputFilename, "generated.bam"])}' \
         '~{bam}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/samtools:v1.7.0_cv3"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.bam"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'SamTools: Sort'
   doc: |-
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/samtools:v1.7.0_cv3

   inputs:
   - id: compression
     label: compression
     doc: |-
       Set the desired compression level for the final output file, ranging from 0 (uncompressed) or 1 (fastest but minimal compression) to 9 (best compression but slowest to write), similarly to gzip(1)'s compression level setting.
       If -l is not used, the default compression level will apply.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -l
   - id: maximumMemory
     label: maximumMemory
     doc: |-
       Approximately the maximum required memory per thread, specified  either in bytes or with a K, M, or G suffix [768 MiB]. To prevent sort from creating a huge number of temporary files, it enforces a minimum value of 1M for this setting.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -m
   - id: sortByReadNames
     label: sortByReadNames
     doc: |-
       Sort by read names (i.e., the QNAME field) rather than by chromosomal coordinates.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -n
   - id: outputType
     label: outputType
     doc: |-
       Write the final output as sam, bam, or cram. By default, samtools tries to select a format based on the -o filename extension; if output is to standard output or no format can be deduced, bam is selected.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -O
   - id: temporaryFilesPrefix
     label: temporaryFilesPrefix
     doc: |-
       Write temporary files to PREFIX.nnnn.bam, or if the specified PREFIX is an existing directory, to PREFIX/samtools.mmm.mmm.tmp.nnnn.bam, where mmm is unique to this invocation of the sort command.
       By default, any temporary files are written alongside the output file, as out.bam.tmp.nnnn.bam, or if output is to standard output, in the current directory as samtools.mmm.mmm.tmp.nnnn.bam.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -T
   - id: threads
     label: threads
     doc: |-
       Set number of sorting and compression threads. By default, operation is single-threaded.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -@
   - id: bam
     label: bam
     type: File
     inputBinding:
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: Output to FILE [stdout].
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: -o
       position: 5

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - samtools
   - sort
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: SamToolsSort


