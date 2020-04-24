:orphan:

GATK4: Base Recalibrator
================================================

*1 contributor · 4 versions*

First pass of the base quality score recalibration. Generates a recalibration table based on various covariates. 
The default covariates are read group, reported quality score, machine cycle, and nucleotide context.

This walker generates tables based on specified covariates. It does a by-locus traversal operating only at sites 
that are in the known sites VCF. ExAc, gnomAD, or dbSNP resources can be used as known sites of variation. 
We assume that all reference mismatches we see are therefore errors and indicative of poor base quality. 
Since there is a large amount of data one can then calculate an empirical probability of error given the 
particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a 
table (of the several covariate values, num observations, num mismatches, empirical quality score).


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.baserecalibrator.versions import Gatk4BaseRecalibrator_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4baserecalibrator_step",
           Gatk4BaseRecalibrator_4_1_3(
               bam=None,
               knownSites=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4baserecalibrator_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4BaseRecalibrator:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4BaseRecalibrator > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       knownSites:
       - knownSites_0.vcf.gz
       - knownSites_1.vcf.gz
       reference: reference.fasta




5. Run Gatk4BaseRecalibrator with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4BaseRecalibrator





Information
------------


:ID: ``Gatk4BaseRecalibrator``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_bqsr_BaseRecalibrator.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     tsv
======  ======  ===============



Additional configuration (inputs)
---------------------------------

==============  ===========================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name            type                         prefix           position  documentation
==============  ===========================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
bam             IndexedBam                   -I                      6  BAM/SAM/CRAM file containing reads
knownSites      Array<CompressedIndexedVCF>  --known-sites          28  **One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis.** This algorithm treats every reference mismatch as an indication of error. However, real genetic variation is expected to mismatch the reference, so it is critical that a database of known polymorphic sites is given to the tool in order to skip over those sites. This tool accepts any number of Feature-containing files (VCF, BCF, BED, etc.) for use as this database. For users wishing to exclude an interval list of known variation simply use -XL my.interval.list to skip over processing those sites. Please note however that the statistics reported by the tool will not accurately reflected those sites skipped by the -XL argument.
reference       FastaWithIndexes             -R                      5  Reference sequence file
tmpDir          Optional<String>             --tmp-dir                  Temp directory to use.
outputFilename  Optional<Filename>           -O                      8  **The output recalibration table filename to create.** After the header, data records occur one per line until the end of the file. The first several items on a line are the values of the individual covariates and will change depending on which covariates were specified at runtime. The last three items are the data- that is, number of observations for this combination of covariates, number of reference mismatches, and the raw empirical quality score calculated by phred-scaling the mismatch rate. Use '/dev/stdout' to print to standard out.
intervals       Optional<bed>                --intervals                -L (BASE) One or more genomic intervals over which to operate
==============  ===========================  =============  ==========  ===============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
