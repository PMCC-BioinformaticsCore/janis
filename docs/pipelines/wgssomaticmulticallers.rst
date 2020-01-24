:orphan:

WGS Somatic (Multi callers)
====================================================

A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 1 contributor · 1 version

:ID: ``WGSSomaticMultiCallers``
:Python: ``janis_pipelines.wgs_somatic.wgssomatic import WGSSomaticMultiCallers``
:Versions: 1.2.0
:Authors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16
:Required inputs:
   - ``normal_inputs: Array<FastqGzPair>``

   - ``tumor_inputs: Array<FastqGzPair>``

   - ``cutadapt_adapters: File``

   - ``gridss_blacklist: bed``

   - ``gatk_intervals: Array<bed>``

   - ``vardict_intervals: Array<bed>``

   - ``vardict_header_lines: File``

   - ``reference: FastaWithDict``

   - ``snps_dbsnp: CompressedIndexedVCF``

   - ``snps_1000gp: CompressedIndexedVCF``

   - ``known_indels: CompressedIndexedVCF``

   - ``mills_indels: CompressedIndexedVCF``
:Outputs: 
   - ``normal_report: Array<Array<Zip>>``

   - ``tumor_report: Array<Array<Zip>>``

   - ``normal_bam: BamPair``

   - ``tumor_bam: BamPair``

   - ``gridss_assembly: VCF``

   - ``variants_gatk: VCF``

   - ``variants_strelka: VCF``

   - ``variants_vardict: VCF``

   - ``variants_gridss: VCF``

   - ``variants_combined: VCF``

Documentation
-------------

URL: *No URL to the documentation was provided*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Embedded Tools
***************

==============================  ======================================================================================================================================
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10d24a2b0>>``
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10d25b128>>``
GATK4 Somatic Variant Caller    ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs              ``Gatk4GatherVcfs/4.1.3.0``
Strelka Somatic Variant Caller  ``strelkaSomaticVariantCaller/v0.1.0``
Gridss                          ``gridss/v2.5.1-dev``
Vardict Somatic Variant Caller  ``vardictSomaticVariantCaller/v0.1.0``
Combine Variants                ``combinevariants/0.0.4``
BCFTools: Sort                  ``bcftoolssort/v1.9``
==============================  ======================================================================================================================================

------

Additional configuration (inputs)
---------------------------------

========================  =======================  ===============
name                      type                     documentation
========================  =======================  ===============
normal_inputs             Array<FastqGzPair>
tumor_inputs              Array<FastqGzPair>
cutadapt_adapters         File
gridss_blacklist          bed
gatk_intervals            Array<bed>
vardict_intervals         Array<bed>
vardict_header_lines      File
reference                 FastaWithDict
snps_dbsnp                CompressedIndexedVCF
snps_1000gp               CompressedIndexedVCF
known_indels              CompressedIndexedVCF
mills_indels              CompressedIndexedVCF
normal_name               Optional<String>
tumor_name                Optional<String>
strelka_intervals         Optional<BedTABIX>
allele_freq_threshold     Optional<Float>
combine_variants_type     Optional<String>
combine_variants_columns  Optional<Array<String>>
========================  =======================  ===============

.
