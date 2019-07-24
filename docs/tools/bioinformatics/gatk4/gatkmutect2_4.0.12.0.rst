:orphan:


GATK4: MuTect2
============================

Description
-------------

Tool identifier: ``gatkmutect2``

Tool path: ``janis_bioinformatics.tools.gatk4.mutect2.latest import Gatk4Mutect2Latest``

Version: 4.0.12.0

Docker: ``broadinstitute/gatk:4.0.12.0``



Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>`_

Tool documentation
******************
Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV) 
and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic 
genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 
for a step-by-step description of the workflow and Article#11127 for an overview of what traditional 
somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory. 
Although we present the tool for somatic calling, it may apply to other contexts, 
such as mitochondrial variant calling.

Outputs
-------
======  ==========  =================
name    type        documentation
======  ==========  =================
out     vcf-gz-tbi  To determine type
======  ==========  =================

Inputs
------
Find the inputs below

Required inputs
***************

==========  =============  ========  ==========  ======================================================================================
name        type           prefix      position  documentation
==========  =============  ========  ==========  ======================================================================================
tumor       BamPair        -I                 6  BAM/SAM/CRAM file containing reads
tumorName   String         -tumor             6  BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode.
normal      BamPair        -I                 5  BAM/SAM/CRAM file containing reads
normalName  String         -normal            6  BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode.
reference   FastaWithDict  -R                 8  Reference sequence file
==========  =============  ========  ==========  ======================================================================================

Optional inputs
***************

========================  ==================  ===============================  ==========  ==============================================================================================================================================================
name                      type                prefix                             position  documentation
========================  ==================  ===============================  ==========  ==============================================================================================================================================================
intervals                 Optional<bed>       -L                                        7  One or more genomic intervals over which to operate
outputFilename            Optional<Filename>  -O                                       20
germlineResource          Optional<VCFIDX>    --germline-resource                      10
afOfAllelesNotInResource  Optional<Float>     --af-of-alleles-not-in-resource          11  Population allele fraction assigned to alleles not found in germline resource. Please see docs/mutect/mutect2.pdf fora derivation of the default value.
panelOfNormals            Optional<VCFIDX>    --panel-of-normals                       10  A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
========================  ==================  ===============================  ==========  ==============================================================================================================================================================


Metadata
********

Author: Michael Franklin


*GATK4: MuTect2 was last updated on 2019-01-24*.
*This page was automatically generated on 2019-07-24*.
