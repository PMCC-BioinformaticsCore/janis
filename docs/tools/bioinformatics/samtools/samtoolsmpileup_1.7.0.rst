:orphan:

SamTools: Mpileup
===================================

*1 contributor · 2 versions*

Generate text pileup output for one or multiple BAM files. Each input file produces a separate group of pileup columns in the output.

Samtools mpileup can still produce VCF and BCF output (with -g or -u), but this feature is deprecated and will be removed in a future release. Please use bcftools mpileup for this instead. (Documentation on the deprecated options has been removed from this manual page, but older versions are available online at <http://www.htslib.org/doc/>.)

Note that there are two orthogonal ways to specify locations in the input file; via -r region and -l file. The former uses (and requires) an index to do random access while the latter streams through the file contents filtering out the specified regions, requiring no index. The two may be used in conjunction. For example a BED file containing locations of genes in chromosome 20 could be specified using -r 20 -l chr20.bed, meaning that the index is used to find chromosome 20 and then it is filtered for the regions listed in the bed file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.samtools.mpileup.versions import SamToolsMpileup_1_7

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "samtoolsmpileup_step",
           SamToolsMpileup_1_7(
               bam=None,
           )
       )
       wf.output("out", source=samtoolsmpileup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for SamToolsMpileup:

.. code-block:: bash

   # user inputs
   janis inputs SamToolsMpileup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam




5. Run SamToolsMpileup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       SamToolsMpileup





Information
------------


:ID: ``SamToolsMpileup``
:URL: `http://www.htslib.org/doc/samtools-mpileup.html <http://www.htslib.org/doc/samtools-mpileup.html>`_
:Versions: 1.9.0, 1.7.0
:Container: biocontainers/samtools:v1.7.0_cv3
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-05-19
:Updated: 2020-05-19



Outputs
-----------

======  ================  ===============
name    type              documentation
======  ================  ===============
out     stdout<TextFile>
======  ================  ===============



Additional configuration (inputs)
---------------------------------

======================  =================  =================  ==========  ========================================================================
name                    type               prefix               position  documentation
======================  =================  =================  ==========  ========================================================================
bam                     IndexedBam                                    10
illuminaEncoding        Optional<Boolean>  --illumina1.3+                 Assume the quality is in the Illumina 1.3+ encoding.
countOrphans            Optional<Boolean>  --count-orphans                do not discard anomalous read pairs
noBAQ                   Optional<Boolean>  --no-BAQ                       disable BAQ (per-Base Alignment Quality)
adjustMQ                Optional<Integer>  --adjust-MQ                    adjust mapping quality; recommended:50, disable:0 [0]
maxDepth                Optional<Integer>  --max-depth                    max per-file depth; avoids excessive memory usage [8000]
redoBAQ                 Optional<Boolean>  --redo-BAQ                     recalculate BAQ on the fly, ignore existing BQs
fastaRef                Optional<File>     --fasta-ref                    skip unlisted positions (chr pos) or regions (BED)
excludeRG               Optional<File>     --exclude-RG                   exclude read groups listed in FILE
positions               Optional<File>     --positions                    skip unlisted positions (chr pos) or regions (BED)
minBQ                   Optional<Integer>  --min-BQ                       Minimum base quality for a base to be considered [13]
minMQ                   Optional<Integer>  --min-MQ                       skip alignments with mapQ smaller than INT [0]
region                  Optional<String>   --region                       region in which pileup is generated
ignoreRG                Optional<Boolean>  --ignore-RG                    ignore RG tags (one BAM = one sample)
inclFlags               Optional<String>   --incl-flags                   required flags: skip reads with mask bits unset []
exclFlags               Optional<String>   --excl-flags                   filter flags: skip reads with mask bits set [UNMAP,SECONDARY,QCFAIL,DUP]
ignoreOverlaps          Optional<Boolean>  --ignore-overlaps              disable read-pair overlap detection
outputBP                Optional<Boolean>  --output-BP                    output base positions on reads
outputMQ                Optional<Boolean>  --output-MQ                    output mapping quality
outputQNAME             Optional<Boolean>  --output-QNAME                 output read names
allPositions            Optional<Boolean>  -a                             output all positions (including zero depth)
absolutelyAllPositions  Optional<Boolean>                                 output absolutely all positions, including unused ref. sequences
reference               Optional<File>     --reference                    Reference sequence FASTA FILE [null]
======================  =================  =================  ==========  ========================================================================
