:orphan:

BEDTools: genomeCoverageBed
=======================================================

``bedtoolsgenomeCoverageBed`` · *1 contributor · 1 version*

bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome. Note: 1. If using BED/GFF/VCF, the input (-i) file must be grouped by chromosome. A simple sort -k 1,1 in.bed > in.sorted.bed will suffice. Also, if using BED/GFF/VCF, one must provide a genome file via the -g argument. 2. If the input is in BAM (-ibam) format, the BAM file must be sorted by position. Using samtools sort aln.bam aln.sorted will suffice.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bedtools.genomecoveragebed.versions import BedToolsGenomeCoverageBed_2_29_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bedtoolsgenomecoveragebed_step",
           BedToolsGenomeCoverageBed_2_29_2(

           )
       )
       wf.output("out", source=bedtoolsgenomecoveragebed_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bedtoolsgenomeCoverageBed:

.. code-block:: bash

   # user inputs
   janis inputs bedtoolsgenomeCoverageBed > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run bedtoolsgenomeCoverageBed with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bedtoolsgenomeCoverageBed





Information
------------

:ID: ``bedtoolsgenomeCoverageBed``
:URL: `https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html <https://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html>`_
:Versions: v2.29.2
:Container: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-04-01
:Updated: 2020-04-01


Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     stdout<TextFile>
======  ================  ===============


Additional configuration (inputs)
---------------------------------

===============  =================  ==========  ==========  ==========================================================================================================================================================================================================================================================================================================================================================================================================
name             type               prefix      position    documentation
===============  =================  ==========  ==========  ==========================================================================================================================================================================================================================================================================================================================================================================================================
depth            Optional<Boolean>  -d                      Report the depth at each genome position (with one-based coordinates). Default behavior is to report a histogram.
depthZero        Optional<Boolean>  -dz                     Report the depth at each genome position (with zero-based coordinates). Reports only non-zero positions. Default behavior is to report a histogram.
BedGraphFormat   Optional<Boolean>  -bg                     Report depth in BedGraph format. For details, see: genome.ucsc.edu/goldenPath/help/bedgraph.html
BedGraphFormata  Optional<Boolean>  -bga                    Report depth in BedGraph format, as above (-bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0  coverage by applying: 'grep -w 0$' to the output.
split            Optional<Boolean>  -split                  Treat 'split' BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR 'N' and 'D' operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).
strand           Optional<String>   -strand                 (STRING): can be + or -. Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).
pairEnd          Optional<Boolean>  -pc                     Calculate coverage of pair-end fragments. Works for BAM files only
fragmentSize     Optional<Boolean>  -fs                     Force to use provided fragment size instead of read length. Works for BAM files only
du               Optional<Boolean>  -du                     Change strand af the mate read (so both reads from the same strand) useful for strand specific. Works for BAM files only
fivePos          Optional<Boolean>  -5                      Calculate coverage of 5' positions (instead of entire interval).
threePos         Optional<Boolean>  -3                      Calculate coverage of 3' positions (instead of entire interval).
max              Optional<Integer>  -max                    Combine all positions with a depth >= max into a single bin in the histogram. Irrelevant for -d and -bedGraph
scale            Optional<Float>    -scale                  Scale the coverage by a constant factor. Each coverage value is multiplied by this factor before being reported. Useful for normalizing coverage by, e.g., reads per million (RPM). Default is 1.0; i.e., unscaled.
trackline        Optional<Boolean>  -trackline              Adds a UCSC/Genome-Browser track line definition in the first line of the output. - See here for more details about track line definition: http://genome.ucsc.edu/goldenPath/help/bedgraph.html - NOTE: When adding a trackline definition, the output BedGraph can be easily uploaded to the Genome Browser as a custom track, BUT CAN NOT be converted into a BigWig file (w/o removing the first line).
trackopts        Optional<String>   -trackopts              Writes additional track line definition parameters in the first line. - Example: -trackopts 'name="My Track" visibility=2 color=255,30,30' Note the use of single-quotes if you have spaces in your parameters.
inputBam         Optional<BAM>      -ibam                   Input bam file. Note: BAM _must_ be sorted by position. A 'samtools sort <BAM>' should suffice.
inputBed         Optional<File>     -iBed                   Input bed file. Must be grouped by chromosome. A simple 'sort -k 1,1 <BED> > <BED>.sorted' will suffice.
inputFile        Optional<File>     -i                      Input file, can be gff/vcf.
genome           Optional<File>     -g                      Genome file. The genome file should tab delimited and structured as follows: <chromName><TAB><chromSize>.
===============  =================  ==========  ==========  ==========================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bedtoolsgenomeCoverageBed {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? depth
       Boolean? depthZero
       Boolean? BedGraphFormat
       Boolean? BedGraphFormata
       Boolean? split
       String? strand
       Boolean? pairEnd
       Boolean? fragmentSize
       Boolean? du
       Boolean? fivePos
       Boolean? threePos
       Int? max
       Float? scale
       Boolean? trackline
       String? trackopts
       File? inputBam
       File? inputBed
       File? inputFile
       File? genome
     }
     command <<<
       set -e
       genomeCoverageBed \
         ~{if (defined(depth) && select_first([depth])) then "-d" else ""} \
         ~{if (defined(depthZero) && select_first([depthZero])) then "-dz" else ""} \
         ~{if (defined(BedGraphFormat) && select_first([BedGraphFormat])) then "-bg" else ""} \
         ~{if (defined(BedGraphFormata) && select_first([BedGraphFormata])) then "-bga" else ""} \
         ~{if (defined(split) && select_first([split])) then "-split" else ""} \
         ~{if defined(strand) then ("-strand '" + strand + "'") else ""} \
         ~{if (defined(pairEnd) && select_first([pairEnd])) then "-pc" else ""} \
         ~{if (defined(fragmentSize) && select_first([fragmentSize])) then "-fs" else ""} \
         ~{if (defined(du) && select_first([du])) then "-du" else ""} \
         ~{if (defined(fivePos) && select_first([fivePos])) then "-5" else ""} \
         ~{if (defined(threePos) && select_first([threePos])) then "-3" else ""} \
         ~{if defined(max) then ("-max " + max) else ''} \
         ~{if defined(scale) then ("-scale " + scale) else ''} \
         ~{if (defined(trackline) && select_first([trackline])) then "-trackline" else ""} \
         ~{if defined(trackopts) then ("-trackopts '" + trackopts + "'") else ""} \
         ~{if defined(inputBam) then ("-ibam '" + inputBam + "'") else ""} \
         ~{if defined(inputBed) then ("-iBed '" + inputBed + "'") else ""} \
         ~{if defined(inputFile) then ("-i '" + inputFile + "'") else ""} \
         ~{if defined(genome) then ("-g '" + genome + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = stdout()
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BEDTools: genomeCoverageBed'
   doc: |-
     bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome. Note: 1. If using BED/GFF/VCF, the input (-i) file must be grouped by chromosome. A simple sort -k 1,1 in.bed > in.sorted.bed will suffice. Also, if using BED/GFF/VCF, one must provide a genome file via the -g argument. 2. If the input is in BAM (-ibam) format, the BAM file must be sorted by position. Using samtools sort aln.bam aln.sorted will suffice.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0

   inputs:
   - id: depth
     label: depth
     doc: |-
       Report the depth at each genome position (with one-based coordinates). Default behavior is to report a histogram.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -d
   - id: depthZero
     label: depthZero
     doc: |-
       Report the depth at each genome position (with zero-based coordinates). Reports only non-zero positions. Default behavior is to report a histogram.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -dz
   - id: BedGraphFormat
     label: BedGraphFormat
     doc: |-
       Report depth in BedGraph format. For details, see: genome.ucsc.edu/goldenPath/help/bedgraph.html
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -bg
   - id: BedGraphFormata
     label: BedGraphFormata
     doc: |-
       Report depth in BedGraph format, as above (-bg). However with this option, regions with zero coverage are also reported. This allows one to quickly extract all regions of a genome with 0  coverage by applying: 'grep -w 0$' to the output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -bga
   - id: split
     label: split
     doc: |-
       Treat 'split' BAM or BED12 entries as distinct BED intervals when computing coverage. For BAM files, this uses the CIGAR 'N' and 'D' operations to infer the blocks for computing coverage. For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds fields (i.e., columns 10,11,12).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -split
   - id: strand
     label: strand
     doc: |-
       (STRING): can be + or -. Calculate coverage of intervals from a specific strand. With BED files, requires at least 6 columns (strand is column 6).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -strand
   - id: pairEnd
     label: pairEnd
     doc: Calculate coverage of pair-end fragments. Works for BAM files only
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -pc
   - id: fragmentSize
     label: fragmentSize
     doc: |-
       Force to use provided fragment size instead of read length. Works for BAM files only
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -fs
   - id: du
     label: du
     doc: |-
       Change strand af the mate read (so both reads from the same strand) useful for strand specific. Works for BAM files only
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -du
   - id: fivePos
     label: fivePos
     doc: Calculate coverage of 5' positions (instead of entire interval).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-5'
   - id: threePos
     label: threePos
     doc: Calculate coverage of 3' positions (instead of entire interval).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-3'
   - id: max
     label: max
     doc: |-
       Combine all positions with a depth >= max into a single bin in the histogram. Irrelevant for -d and -bedGraph
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -max
   - id: scale
     label: scale
     doc: |-
       Scale the coverage by a constant factor. Each coverage value is multiplied by this factor before being reported. Useful for normalizing coverage by, e.g., reads per million (RPM). Default is 1.0; i.e., unscaled.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -scale
   - id: trackline
     label: trackline
     doc: |-
       Adds a UCSC/Genome-Browser track line definition in the first line of the output. - See here for more details about track line definition: http://genome.ucsc.edu/goldenPath/help/bedgraph.html - NOTE: When adding a trackline definition, the output BedGraph can be easily uploaded to the Genome Browser as a custom track, BUT CAN NOT be converted into a BigWig file (w/o removing the first line).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -trackline
   - id: trackopts
     label: trackopts
     doc: |-
       Writes additional track line definition parameters in the first line. - Example: -trackopts 'name="My Track" visibility=2 color=255,30,30' Note the use of single-quotes if you have spaces in your parameters.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -trackopts
   - id: inputBam
     label: inputBam
     doc: |-
       Input bam file. Note: BAM _must_ be sorted by position. A 'samtools sort <BAM>' should suffice.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -ibam
   - id: inputBed
     label: inputBed
     doc: |-
       Input bed file. Must be grouped by chromosome. A simple 'sort -k 1,1 <BED> > <BED>.sorted' will suffice.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -iBed
   - id: inputFile
     label: inputFile
     doc: Input file, can be gff/vcf.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -i
   - id: genome
     label: genome
     doc: |-
       Genome file. The genome file should tab delimited and structured as follows: <chromName><TAB><chromSize>.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -g

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - genomeCoverageBed
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bedtoolsgenomeCoverageBed


