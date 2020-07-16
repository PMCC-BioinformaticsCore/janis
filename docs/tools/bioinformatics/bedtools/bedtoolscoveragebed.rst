:orphan:

BEDTools: coverageBed
===========================================

*0 contributors · 1 version*

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
:Authors: 
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
