:orphan:

WGS Somatic (GATK only)
========================================

A somatic tumor-normal variant-calling WGS pipeline using only GATK Mutect2 · 1 contributor · 1 version

:ID: ``WGSSomaticGATK``
:Python: ``janis_pipelines.wgs_somatic_gatk.wgssomaticgatk import WGSSomaticGATK``
:Versions: 1.1.0
:Authors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16
:Required inputs:
   - ``normalInputs: Array<FastqGzPair>``

   - ``tumorInputs: Array<FastqGzPair>``

   - ``gatkIntervals: Array<bed>``

   - ``cutadapt_adapters: File``

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

   - ``variants_gatk: CompressedVCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

============================  ======================================================================================================================================
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10d597128>>``
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10d5a0198>>``
GATK4 Somatic Variant Caller  ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs            ``Gatk4GatherVcfs/4.1.3.0``
BCFTools: Sort                ``bcftoolssort/v1.9``
============================  ======================================================================================================================================

------

Additional configuration (inputs)
---------------------------------

=================  ====================  ===============
name               type                  documentation
=================  ====================  ===============
normalInputs       Array<FastqGzPair>
tumorInputs        Array<FastqGzPair>
gatkIntervals      Array<bed>
cutadapt_adapters  File
reference          FastaWithDict
snps_dbsnp         CompressedIndexedVCF
snps_1000gp        CompressedIndexedVCF
known_indels       CompressedIndexedVCF
mills_indels       CompressedIndexedVCF
normalName         Optional<String>
tumorName          Optional<String>
=================  ====================  ===============

.
