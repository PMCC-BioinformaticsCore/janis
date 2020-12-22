:orphan:

Tabix
=============

``tabix`` · *1 contributor · 2 versions*

tabix – Generic indexer for TAB-delimited genome position files

Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
and compressed by bgzip which has a gzip(1) like interface.

After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
"chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

Fast data retrieval also works over network if URI is given as a file name and in this case the 
index file will be downloaded if it is not present locally.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.htslib.tabix.tabix_1_2_1 import Tabix_1_2_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "tabix_step",
           Tabix_1_2_1(
               inp=None,
           )
       )
       wf.output("out", source=tabix_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for tabix:

.. code-block:: bash

   # user inputs
   janis inputs tabix > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp: inp.gz




5. Run tabix with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       tabix





Information
------------

:ID: ``tabix``
:URL: `http://www.htslib.org/doc/tabix.html <http://www.htslib.org/doc/tabix.html>`_
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
out     Gzipped<File>
======  =============  ===============


Additional configuration (inputs)
---------------------------------

===========  =================  ==============  ==========  ==============================================================================================================================================================================================================================================================================================================
name         type               prefix            position  documentation
===========  =================  ==============  ==========  ==============================================================================================================================================================================================================================================================================================================
inp          Gzipped<File>                               8  File from which to create the index. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
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

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task tabix {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File inp
       String? preset
       Boolean? zeroBased
       Int? begin
       String? comment
       Boolean? csi
       Int? end
       Boolean? force
       Int? minShift
       Int? sequence
       Int? skipLines
       Boolean? printHeader
       Boolean? onlyHeader
       Boolean? listChroms
       File? reheader
       File? regions
       File? targets
     }
     command <<<
       set -e
       cp -f '~{inp}' '.'
       tabix \
         ~{if (defined(zeroBased) && select_first([zeroBased])) then "--zero-based" else ""} \
         ~{if (defined(csi) && select_first([csi])) then "--csi" else ""} \
         ~{if (defined(force) && select_first([force])) then "--force" else ""} \
         ~{if defined(minShift) then ("--min-shift " + minShift) else ''} \
         ~{if (defined(printHeader) && select_first([printHeader])) then "--print-header" else ""} \
         ~{if (defined(onlyHeader) && select_first([onlyHeader])) then "--only-header" else ""} \
         ~{if (defined(listChroms) && select_first([listChroms])) then "--list-chroms" else ""} \
         ~{if defined(reheader) then ("--reheader '" + reheader + "'") else ""} \
         ~{if defined(select_first([preset, "vcf"])) then ("--preset '" + select_first([preset, "vcf"]) + "'") else ""} \
         ~{if defined(sequence) then ("--sequence " + sequence) else ''} \
         ~{if defined(begin) then ("--begin " + begin) else ''} \
         ~{if defined(end) then ("--end " + end) else ''} \
         ~{if defined(skipLines) then ("--skip-lines " + skipLines) else ''} \
         ~{if defined(comment) then ("--comment '" + comment + "'") else ""} \
         '~{basename(inp)}' \
         ~{if defined(regions) then ("--regions '" + regions + "'") else ""} \
         ~{if defined(targets) then ("--targets '" + targets + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biodckrdev/htslib:1.2.1"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = basename(inp)
       File out_tbi = basename(inp) + ".tbi"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Tabix
   doc: |-
     tabix – Generic indexer for TAB-delimited genome position files

     Tabix indexes a TAB-delimited genome position file in.tab.bgz and creates an index file (in.tab.bgz.tbi or 
     in.tab.bgz.csi) when region is absent from the command-line. The input data file must be position sorted 
     and compressed by bgzip which has a gzip(1) like interface.

     After indexing, tabix is able to quickly retrieve data lines overlapping regions specified in the format 
     "chr:beginPos-endPos". (Coordinates specified in this region format are 1-based and inclusive.)

     Fast data retrieval also works over network if URI is given as a file name and in this case the 
     index file will be downloaded if it is not present locally.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: InitialWorkDirRequirement
     listing:
     - entry: $(inputs.inp)
   - class: DockerRequirement
     dockerPull: biodckrdev/htslib:1.2.1

   inputs:
   - id: inp
     label: inp
     doc: |-
       File from which to create the index. The input data file must be position sorted and compressed by bgzip which has a gzip(1) like interface.
     type: File
     inputBinding:
       position: 8
   - id: preset
     label: preset
     doc: |-
       -p: Input format for indexing. Valid values are: gff, bed, sam, vcf. This option should not be applied together with any of -s, -b, -e, -c and -0; it is not used for data retrieval because this setting is stored in the index file. [gff]
     type: string
     default: vcf
     inputBinding:
       prefix: --preset
       position: 2
   - id: zeroBased
     label: zeroBased
     doc: |-
       -0: Specify that the position in the data file is 0-based (e.g. UCSC files) rather than 1-based.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --zero-based
       position: 1
   - id: begin
     label: begin
     doc: '-b: Column of start chromosomal position. [4]'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --begin
       position: 4
   - id: comment
     label: comment
     doc: '-c: Skip lines started with character CHAR. [#]'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --comment
       position: 7
   - id: csi
     label: csi
     doc: '-C: Produce CSI format index instead of classical tabix or BAI style indices.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --csi
       position: 1
   - id: end
     label: end
     doc: |-
       -e: Column of end chromosomal position. The end column can be the same as the start column. [5]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --end
       position: 5
   - id: force
     label: force
     doc: '-f: Force to overwrite the index file if it is present.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force
       position: 1
   - id: minShift
     label: minShift
     doc: '-m: set minimal interval size for CSI indices to 2^INT [14]'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-shift
       position: 1
   - id: sequence
     label: sequence
     doc: |-
       -s: Column of sequence name. Option -s, -b, -e, -S, -c and -0 are all stored in the index file and thus not used in data retrieval. [1]
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --sequence
       position: 3
   - id: skipLines
     label: skipLines
     doc: '-S: Skip first INT lines in the data file. [0]'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --skip-lines
       position: 6
   - id: printHeader
     label: printHeader
     doc: '-h: Print also the header/meta lines.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --print-header
       position: 1
   - id: onlyHeader
     label: onlyHeader
     doc: '-H: Print only the header/meta lines.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --only-header
       position: 1
   - id: listChroms
     label: listChroms
     doc: '-l: List the sequence names stored in the index file.'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --list-chroms
       position: 1
   - id: reheader
     label: reheader
     doc: '-r: Replace the header with the content of FILE'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --reheader
       position: 1
   - id: regions
     label: regions
     doc: |-
       -R: Restrict to regions listed in the FILE. The FILE can be BED file (requires .bed, .bed.gz, .bed.bgz file name extension) or a TAB-delimited file with CHROM, POS, and, optionally, POS_TO columns, where positions are 1-based and inclusive. When this option is in use, the input file may not be sorted.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --regions
       position: 11
   - id: targets
     label: targets
     doc: |-
       -T: Similar to -R but the entire input will be read sequentially and regions not listed in FILE will be skipped
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --targets
       position: 11

   outputs:
   - id: out
     label: out
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $(inputs.inp.basename)
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: tabix
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: tabix


