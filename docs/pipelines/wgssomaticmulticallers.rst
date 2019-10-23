:orphan:

WGS Somatic (Multi callers)
====================================================

A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 1 contributor · 1 version

:ID: ``WGSSomaticMultiCallers``
:Python: ``janis_pipelines.wgs_somatic.wgssomatic import WGSSomaticMultiCallers``
:Versions: v0.1.0
:Contributors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16
:Required inputs:
   - ``normalInputs: Array<FastqGzPair>``

   - ``tumorInputs: Array<FastqGzPair>``

   - ``gatkIntervals: Array<bed>``

   - ``vardictIntervals: Array<bed>``

   - ``vardictHeaderLines: File``

   - ``reference: FastaWithDict``

   - ``snps_dbsnp: CompressedIndexedVCF``

   - ``snps_1000gp: CompressedIndexedVCF``

   - ``known_indels: CompressedIndexedVCF``

   - ``mills_indels: CompressedIndexedVCF``
:Outputs: 
   - ``normalBam: BamPair``

   - ``tumorBam: BamPair``

   - ``normalReport: Array<Array<Zip>>``

   - ``tumorReport: Array<Array<Zip>>``

   - ``variants_gatk: VCF``

   - ``variants_strelka: VCF``

   - ``variants_vardict: VCF``

   - ``variants_combined: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

==============================  ======================================
                                ``somatic_subpipeline/v0.1.0``
GATK4 Somatic Variant Caller    ``GATK4_SomaticVariantCaller/v0.1.0``
GATK4: Gather VCFs              ``Gatk4GatherVcfs/4.1.3.0``
Strelka Somatic Variant Caller  ``strelkaSomaticVariantCaller/v0.1.0``
Vardict Somatic Variant Caller  ``vardictSomaticVariantCaller/v0.1.0``
Combine Variants                ``combinevariants/0.0.4``
BCFTools: Sort                  ``bcftoolssort/v1.9``
==============================  ======================================

------

Inputs
------

=======================  =======================  ===============
name                     type                     documentation
=======================  =======================  ===============
normalInputs             Array<FastqGzPair>
tumorInputs              Array<FastqGzPair>
gatkIntervals            Array<bed>
vardictIntervals         Array<bed>
vardictHeaderLines       File
reference                FastaWithDict
snps_dbsnp               CompressedIndexedVCF
snps_1000gp              CompressedIndexedVCF
known_indels             CompressedIndexedVCF
mills_indels             CompressedIndexedVCF
normalName               Optional<String>
tumorName                Optional<String>
strelkaIntervals         Optional<BedTABIX>
alleleFreqThreshold      Optional<Float>
combineVariants_type     Optional<String>
combineVariants_columns  Optional<Array<String>>
=======================  =======================  ===============

.
