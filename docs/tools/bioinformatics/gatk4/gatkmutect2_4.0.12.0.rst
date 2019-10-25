:orphan:

GATK4: MuTect2
============================

1 contributor · 4 versions

:ID: ``gatkmutect2``
:Python: ``janis_bioinformatics.tools.gatk4.mutect2.versions import GatkMutect2_4_0``
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.0.12.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24
:Required inputs:
   - ``tumor: BamPair``

   - ``tumorName: String``

   - ``normal: BamPair``

   - ``normalName: String``

   - ``reference: FastaWithDict``
:Outputs: 
   - ``out: CompressedIndexedVCF``

Documentation
-------------

URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.10.0/org_broadinstitute_hellbender_tools_walkers_mutect_Mutect2.php>`_

Call somatic short variants via local assembly of haplotypes. Short variants include single nucleotide (SNV) 
and insertion and deletion (indel) variants. The caller combines the DREAM challenge-winning somatic 
genotyping engine of the original MuTect (Cibulskis et al., 2013) with the assembly-based machinery of HaplotypeCaller.

This tool is featured in the Somatic Short Mutation calling Best Practice Workflow. See Tutorial#11136 
for a step-by-step description of the workflow and Article#11127 for an overview of what traditional 
somatic calling entails. For the latest pipeline scripts, see the Mutect2 WDL scripts directory. 
Although we present the tool for somatic calling, it may apply to other contexts, 
such as mitochondrial variant calling.

------

Additional configuration (inputs)
---------------------------------

========================  ====================  ==============================================================================================================================================================
name                      type                  documentation
========================  ====================  ==============================================================================================================================================================
tumor                     BamPair               BAM/SAM/CRAM file containing reads
tumorName                 String                BAM sample name of tumor. May be URL-encoded as output by GetSampleName with -encode.
normal                    BamPair               BAM/SAM/CRAM file containing reads
normalName                String                BAM sample name of normal. May be URL-encoded as output by GetSampleName with -encode.
reference                 FastaWithDict         Reference sequence file
intervals                 Optional<bed>         One or more genomic intervals over which to operate
outputFilename            Optional<Filename>
germlineResource          Optional<IndexedVCF>
afOfAllelesNotInResource  Optional<Float>       Population allele fraction assigned to alleles not found in germline resource. Please see docs/mutect/mutect2.pdf fora derivation of the default value.
panelOfNormals            Optional<IndexedVCF>  A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
========================  ====================  ==============================================================================================================================================================

