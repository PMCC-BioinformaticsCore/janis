:orphan:

GATK4: GetPileupSummaries
===================================================

1 contributor · 3 versions

:ID: ``GATK4GetPileupSummaries``
:Python: ``janis_bioinformatics.tools.gatk4.getpileupsummaries.versions import Gatk4GetPileUpSummaries_4_1_4``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09
:Required inputs:
   - ``bam: Array<BamPair>``

   - ``sites: CompressedIndexedVCF``
:Outputs: 
   - ``out: TextFile``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/org_broadinstitute_hellbender_tools_walkers_contamination_GetPileupSummaries.php>`_

Summarizes counts of reads that support reference, alternate and other alleles for given sites. Results can be used with CalculateContamination.
The tool requires a common germline variant sites VCF, e.g. the gnomAD resource, with population allele frequencies (AF) in the INFO field. This resource must contain only biallelic SNPs and can be an eight-column sites-only VCF. The tool ignores the filter status of the sites. See the GATK Resource Bundle for an example human file.

------

Additional configuration (inputs)
---------------------------------

==============  ==============================  =============================================================
name            type                            documentation
==============  ==============================  =============================================================
bam             Array<BamPair>                  The SAM/BAM/CRAM file containing reads.
sites           CompressedIndexedVCF            sites of common biallelic variants
intervals       Optional<CompressedIndexedVCF>  -L (BASE) One or more genomic intervals over which to operate
pileupTableOut  Optional<Filename>
==============  ==============================  =============================================================

