:orphan:

BEDTools: coverageBed
===========================================

``bedtoolsCoverageBed`` · *1 contributor · 1 version*

The bedtools coverage tool computes both the depth and breadth of coverage of features in file B on the features in file A. For example, bedtools coverage can compute the coverage of sequence alignments (file B) across 1 kilobase (arbitrary) windows (file A) tiling a genome of interest. One advantage that bedtools coverage offers is that it not only counts the number of features that overlap an interval in file A, it also computes the fraction of bases in the interval in A that were overlapped by one or more features. Thus, bedtools coverage also computes the breadth of coverage observed for each interval in A.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bedtools.coveragebed.versions import BedToolsCoverageBed_2_29_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bedtoolscoveragebed_step",
           BedToolsCoverageBed_2_29_2(
               inputABed=None,
               inputBBam=None,
           )
       )
       wf.output("out", source=bedtoolscoveragebed_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bedtoolsCoverageBed:

.. code-block:: bash

   # user inputs
   janis inputs bedtoolsCoverageBed > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputABed: inputABed.bed
       inputBBam: inputBBam.bam




5. Run bedtoolsCoverageBed with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bedtoolsCoverageBed





Information
------------

:ID: ``bedtoolsCoverageBed``
:URL: `https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html <https://bedtools.readthedocs.io/en/latest/content/tools/coverage.html>`_
:Versions: v2.29.2
:Container: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-02-20
:Updated: 2020-02-26


Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     stdout<TextFile>
======  ================  ===============


Additional configuration (inputs)
---------------------------------

=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================
name                   type               prefix        position    documentation
=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================
inputABed              bed                -a                        input file a: only bed is supported. May be followed with multiple databases and/or  wildcard (*) character(s).
inputBBam              BAM                -b                        input file b: only bam is supported.
strandedness           Optional<Boolean>  -s                        Require same strandedness.  That is, only report hits in B that overlap A on the _same_ strand. - By default, overlaps are reported without respect to strand.
differentStrandedness  Optional<Boolean>  -S                        Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. - By default, overlaps are reported without respect to strand.
fractionA              Optional<Float>    -f                        Minimum overlap required as a fraction of A. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
fractionB              Optional<Float>    -F                        Minimum overlap required as a fraction of B. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
reciprocalFraction     Optional<Boolean>  -r                        Require that the fraction overlap be reciprocal for A AND B. - In other words, if -f is 0.90 and -r is used, this requires that B overlap 90% of A and A _also_ overlaps 90% of B.
minFraction            Optional<Boolean>  -r                        Require that the minimum fraction be satisfied for A OR B. - In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered. Without -e, both fractions would have to be satisfied.
split                  Optional<Boolean>  -split                    Treat 'split' BAM or BED12 entries as distinct BED intervals.
genome                 Optional<File>     -g                        Provide a genome file to enforce consistent chromosome sort order across input files. Only applies when used with -sorted option.
noNameCheck            Optional<Boolean>  -nonamecheck              For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. 'chr1' vs 'chr01'.
sorted                 Optional<Boolean>  -sorted                   Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
header                 Optional<Boolean>  -header                   Print the header from the A file prior to results.
noBuf                  Optional<Boolean>  -nobuf                    Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
bufMem                 Optional<Integer>  -iobuf                    Specify amount of memory to use for input buffer. Takes an integer argument. Optional suffixes K/M/G supported. Note: currently has no effect with compressed files.
histogram              Optional<Boolean>  -hist                     Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A. Output (tab delimited) after each feature in A: 1) depth 2) # bases at depth 3) size of A 4) % of A at depth.
depth                  Optional<Boolean>  -d                        Report the depth at each position in each A feature. Positions reported are one based.  Each position and depth follow the complete A feature.
counts                 Optional<Boolean>  -counts                   Only report the count of overlaps, don't compute fraction, etc.
mean                   Optional<Boolean>  -mean                     Report the mean depth of all positions in each A feature.
=====================  =================  ============  ==========  ===========================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bedtoolsCoverageBed {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Boolean? strandedness
       Boolean? differentStrandedness
       Float? fractionA
       Float? fractionB
       Boolean? reciprocalFraction
       Boolean? minFraction
       Boolean? split
       File? genome
       Boolean? noNameCheck
       Boolean? sorted
       Boolean? header
       Boolean? noBuf
       Int? bufMem
       File inputABed
       File inputBBam
       Boolean? histogram
       Boolean? depth
       Boolean? counts
       Boolean? mean
     }
     command <<<
       set -e
       coverageBed \
         ~{if (defined(strandedness) && select_first([strandedness])) then "-s" else ""} \
         ~{if (defined(differentStrandedness) && select_first([differentStrandedness])) then "-S" else ""} \
         ~{if defined(fractionA) then ("-f " + fractionA) else ''} \
         ~{if defined(fractionB) then ("-F " + fractionB) else ''} \
         ~{if (defined(reciprocalFraction) && select_first([reciprocalFraction])) then "-r" else ""} \
         ~{if (defined(minFraction) && select_first([minFraction])) then "-r" else ""} \
         ~{if (defined(split) && select_first([split])) then "-split" else ""} \
         ~{if defined(genome) then ("-g '" + genome + "'") else ""} \
         ~{if (defined(noNameCheck) && select_first([noNameCheck])) then "-nonamecheck" else ""} \
         ~{if (defined(sorted) && select_first([sorted])) then "-sorted" else ""} \
         ~{if (defined(header) && select_first([header])) then "-header" else ""} \
         ~{if (defined(noBuf) && select_first([noBuf])) then "-nobuf" else ""} \
         ~{if defined(bufMem) then ("-iobuf " + bufMem) else ''} \
         -a '~{inputABed}' \
         -b '~{inputBBam}' \
         ~{if (defined(histogram) && select_first([histogram])) then "-hist" else ""} \
         ~{if (defined(depth) && select_first([depth])) then "-d" else ""} \
         ~{if (defined(counts) && select_first([counts])) then "-counts" else ""} \
         ~{if (defined(mean) && select_first([mean])) then "-mean" else ""}
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
   label: 'BEDTools: coverageBed'
   doc: |-
     The bedtools coverage tool computes both the depth and breadth of coverage of features in file B on the features in file A. For example, bedtools coverage can compute the coverage of sequence alignments (file B) across 1 kilobase (arbitrary) windows (file A) tiling a genome of interest. One advantage that bedtools coverage offers is that it not only counts the number of features that overlap an interval in file A, it also computes the fraction of bases in the interval in A that were overlapped by one or more features. Thus, bedtools coverage also computes the breadth of coverage observed for each interval in A.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0

   inputs:
   - id: strandedness
     label: strandedness
     doc: |-
       Require same strandedness.  That is, only report hits in B that overlap A on the _same_ strand. - By default, overlaps are reported without respect to strand.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -s
   - id: differentStrandedness
     label: differentStrandedness
     doc: |-
       Require different strandedness.  That is, only report hits in B that overlap A on the _opposite_ strand. - By default, overlaps are reported without respect to strand.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -S
   - id: fractionA
     label: fractionA
     doc: |-
       Minimum overlap required as a fraction of A. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -f
   - id: fractionB
     label: fractionB
     doc: |-
       Minimum overlap required as a fraction of B. - Default is 1E-9 (i.e., 1bp). - FLOAT (e.g. 0.50)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -F
   - id: reciprocalFraction
     label: reciprocalFraction
     doc: |-
       Require that the fraction overlap be reciprocal for A AND B. - In other words, if -f is 0.90 and -r is used, this requires that B overlap 90% of A and A _also_ overlaps 90% of B.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -r
   - id: minFraction
     label: minFraction
     doc: |-
       Require that the minimum fraction be satisfied for A OR B. - In other words, if -e is used with -f 0.90 and -F 0.10 this requires that either 90% of A is covered OR 10% of  B is covered. Without -e, both fractions would have to be satisfied.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -r
   - id: split
     label: split
     doc: Treat 'split' BAM or BED12 entries as distinct BED intervals.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -split
   - id: genome
     label: genome
     doc: |-
       Provide a genome file to enforce consistent chromosome sort order across input files. Only applies when used with -sorted option.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -g
   - id: noNameCheck
     label: noNameCheck
     doc: |-
       For sorted data, don't throw an error if the file has different naming conventions for the same chromosome. ex. 'chr1' vs 'chr01'.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -nonamecheck
   - id: sorted
     label: sorted
     doc: Use the 'chromsweep' algorithm for sorted (-k1,1 -k2,2n) input.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -sorted
   - id: header
     label: header
     doc: Print the header from the A file prior to results.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -header
   - id: noBuf
     label: noBuf
     doc: |-
       Disable buffered output. Using this option will cause each line of output to be printed as it is generated, rather than saved in a buffer. This will make printing large output files noticeably slower, but can be useful in conjunction with other software tools and scripts that need to process one line of bedtools output at a time.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -nobuf
   - id: bufMem
     label: bufMem
     doc: |-
       Specify amount of memory to use for input buffer. Takes an integer argument. Optional suffixes K/M/G supported. Note: currently has no effect with compressed files.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -iobuf
   - id: inputABed
     label: inputABed
     doc: |-
       input file a: only bed is supported. May be followed with multiple databases and/or  wildcard (*) character(s). 
     type: File
     inputBinding:
       prefix: -a
   - id: inputBBam
     label: inputBBam
     doc: 'input file b: only bam is supported.'
     type: File
     inputBinding:
       prefix: -b
   - id: histogram
     label: histogram
     doc: |-
       Report a histogram of coverage for each feature in A as well as a summary histogram for _all_ features in A. Output (tab delimited) after each feature in A: 1) depth 2) # bases at depth 3) size of A 4) % of A at depth.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -hist
   - id: depth
     label: depth
     doc: |-
       Report the depth at each position in each A feature. Positions reported are one based.  Each position and depth follow the complete A feature.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -d
   - id: counts
     label: counts
     doc: Only report the count of overlaps, don't compute fraction, etc.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -counts
   - id: mean
     label: mean
     doc: Report the mean depth of all positions in each A feature.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -mean

   outputs:
   - id: out
     label: out
     type: stdout
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - coverageBed
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bedtoolsCoverageBed


