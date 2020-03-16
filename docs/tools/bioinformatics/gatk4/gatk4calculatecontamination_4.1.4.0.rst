:orphan:

GATK4: CalculateContamination
===========================================================

*1 contributor · 3 versions*

:ID: ``Gatk4CalculateContamination``
:Python: ``janis_bioinformatics.tools.gatk4.calculatecontaminations.versions import Gatk4CalculateContamination_4_1_4``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0
:Container: broadinstitute/gatk:4.1.4.0
:Authors: Hollizeck Sebastian
:Citations: TBD
:Created: 2019-09-09
:Updated: 2019-09-09
:Required inputs:
   - ``pileupTable: File``
:Outputs: 
   - ``contOut: File``

   - ``segOut: File``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.2.0/org_broadinstitute_hellbender_tools_walkers_contamination_CalculateContamination.php>`_

Calculates the fraction of reads coming from cross-sample contamination, given results from GetPileupSummaries. The resulting contamination table is used with FilterMutectCalls.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 for a step-by-step description of the workflow and Article#11127 for an overview of what traditional somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory.

This tool borrows from ContEst by Cibulskis et al the idea of estimating contamination from ref reads at hom alt sites. However, ContEst uses a probabilistic model that assumes a diploid genotype with no copy number variation and independent contaminating reads. That is, ContEst assumes that each contaminating read is drawn randomly and independently from a different human. This tool uses a simpler estimate of contamination that relaxes these assumptions. In particular, it works in the presence of copy number variations and with an arbitrary number of contaminating samples. In addition, this tool is designed to work well with no matched normal data. However, one can run GetPileupSummaries on a matched normal bam file and input the result to this tool.

------

None

Additional configuration (inputs)
---------------------------------

====================  ==================  ==================================  ==========  =============================================================================================================================================
name                  type                prefix                                position  documentation
====================  ==================  ==================================  ==========  =============================================================================================================================================
pileupTable           File                -I                                              pileup table from summarize pileup
contaminationTable    Optional<File>      --contamination-table                           Tables containing contamination information.
segmentationFile      Optional<File>      --tumor-segmentation                            Tables containing tumor segments' minor allele fractions for germline hets emitted by CalculateContamination
statsFile             Optional<File>      --stats                                         The Mutect stats file output by Mutect2
readOrientationModel  Optional<File>      --orientation-bias-artifact-priors              One or more .tar.gz files containing tables of prior artifact probabilities for the read orientation filter model, one table per tumor sample
segmentationFileOut   Optional<Filename>  --tumor-segmentation                            Reference sequence file
contaminationFileOut  Optional<Filename>  -O                                           2
====================  ==================  ==================================  ==========  =============================================================================================================================================

