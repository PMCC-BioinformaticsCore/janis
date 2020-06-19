:orphan:

GATK4: GetPileupSummaries
========================================================

*1 contributor · 3 versions*

Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.getpileupsummaries.versions import Gatk4GetPileUpSummariesCram_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4getpileupsummaries_cram_step",
           Gatk4GetPileUpSummariesCram_4_1_3(
               bam=None,
               sites=None,
           )
       )
       wf.output("out", source=gatk4getpileupsummaries_cram_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GetPileupSummaries_cram:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GetPileupSummaries_cram > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam:
       - bam_0.cram
       - bam_1.cram
       sites: sites.vcf.gz




5. Run Gatk4GetPileupSummaries_cram with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GetPileupSummaries_cram





Information
------------


:ID: ``Gatk4GetPileupSummaries_cram``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.3.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09



Outputs
-----------

======  ========  ================================
name    type      documentation
======  ========  ================================
out     TextFile  Table containing the pileup info
======  ========  ================================



Additional configuration (inputs)
---------------------------------

==============  ==============================  ===========  ==========  =============================================================
name            type                            prefix         position  documentation
==============  ==============================  ===========  ==========  =============================================================
bam             Array<CramPair>                 -I                    0  The SAM/BAM/CRAM file containing reads.
sites           CompressedIndexedVCF            -V                       sites of common biallelic variants
intervals       Optional<CompressedIndexedVCF>  --intervals              -L (BASE) One or more genomic intervals over which to operate
pileupTableOut  Optional<Filename>              -O                    1
==============  ==============================  ===========  ==========  =============================================================
