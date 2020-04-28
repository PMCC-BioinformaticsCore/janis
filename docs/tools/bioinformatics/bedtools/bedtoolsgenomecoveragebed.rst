:orphan:

BEDTools: genomeCoverageBed
=======================================================

*0 contributors · 1 version*

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
:Authors: 
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
