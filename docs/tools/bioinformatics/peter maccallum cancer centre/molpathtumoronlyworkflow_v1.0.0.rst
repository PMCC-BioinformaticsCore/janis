:orphan:

Molpath Tumor Only Workflow
======================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.molpathTumorOnlyWorkflow import MolpathTumorOnly_1_0_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "molpathtumoronlyworkflow_step",
           MolpathTumorOnly_1_0_0(
               sample_name=None,
               fastqs=None,
               seqrun=None,
               reference=None,
               region_bed=None,
               region_bed_extended=None,
               region_bed_annotated=None,
               genecoverage_bed=None,
               genome_file=None,
               panel_name=None,
               vcfcols=None,
               snps_dbsnp=None,
               snps_1000gp=None,
               known_indels=None,
               mills_indels=None,
               mutalyzer_server=None,
               pathos_db=None,
               maxRecordsInRam=None,
               gnomad=None,
           )
       )
       wf.output("fastq_qc", source=molpathtumoronlyworkflow_step.fastq_qc)
       wf.output("markdups_bam", source=molpathtumoronlyworkflow_step.markdups_bam)
       wf.output("doc_out", source=molpathtumoronlyworkflow_step.doc_out)
       wf.output("summary", source=molpathtumoronlyworkflow_step.summary)
       wf.output("gene_summary", source=molpathtumoronlyworkflow_step.gene_summary)
       wf.output("region_summary", source=molpathtumoronlyworkflow_step.region_summary)
       wf.output("gridss_vcf", source=molpathtumoronlyworkflow_step.gridss_vcf)
       wf.output("gridss_bam", source=molpathtumoronlyworkflow_step.gridss_bam)
       wf.output("haplotypecaller_vcf", source=molpathtumoronlyworkflow_step.haplotypecaller_vcf)
       wf.output("haplotypecaller_bam", source=molpathtumoronlyworkflow_step.haplotypecaller_bam)
       wf.output("haplotypecaller_norm", source=molpathtumoronlyworkflow_step.haplotypecaller_norm)
       wf.output("mutect2_vcf", source=molpathtumoronlyworkflow_step.mutect2_vcf)
       wf.output("mutect2_bam", source=molpathtumoronlyworkflow_step.mutect2_bam)
       wf.output("mutect2_norm", source=molpathtumoronlyworkflow_step.mutect2_norm)
       wf.output("addbamstats_vcf", source=molpathtumoronlyworkflow_step.addbamstats_vcf)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for MolpathTumorOnlyWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs MolpathTumorOnlyWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastqs:
       - - fastqs_0.fastq.gz
         - fastqs_1.fastq.gz
       - - fastqs_0.fastq.gz
         - fastqs_1.fastq.gz
       genecoverage_bed: genecoverage_bed.bed
       genome_file: genome_file.txt
       gnomad: gnomad.vcf.gz
       known_indels: known_indels.vcf.gz
       maxRecordsInRam: 0
       mills_indels: mills_indels.vcf.gz
       mutalyzer_server: <value>
       panel_name: <value>
       pathos_db: <value>
       reference: reference.fasta
       region_bed: region_bed.bed
       region_bed_annotated: region_bed_annotated.bed
       region_bed_extended: region_bed_extended.bed
       sample_name: <value>
       seqrun: <value>
       snps_1000gp: snps_1000gp.vcf.gz
       snps_dbsnp: snps_dbsnp.vcf.gz
       vcfcols: vcfcols.txt




5. Run MolpathTumorOnlyWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       MolpathTumorOnlyWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``MolpathTumorOnlyWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

====================  ====================  ===============
name                  type                  documentation
====================  ====================  ===============
fastq_qc              Array<Array<Zip>>
markdups_bam          IndexedBam
doc_out               TextFile
summary               csv
gene_summary          TextFile
region_summary        TextFile
gridss_vcf            VCF
gridss_bam            BAM
haplotypecaller_vcf   CompressedIndexedVCF
haplotypecaller_bam   IndexedBam
haplotypecaller_norm  VCF
mutect2_vcf           CompressedIndexedVCF
mutect2_bam           IndexedBam
mutect2_norm          VCF
addbamstats_vcf       VCF
====================  ====================  ===============


Embedded Tools
***************

======================================================================  ======================================================
FastQC                                                                  ``fastqc/v0.11.5``
Parse FastQC Adaptors                                                   ``ParseFastqcAdaptors/v0.1.0``
Align and sort reads                                                    ``BwaAligner/1.0.0``
Merge and Mark Duplicates                                               ``mergeAndMarkBams/4.1.3``
Annotate GATK3 DepthOfCoverage Workflow                                 ``AnnotateDepthOfCoverage/v0.1.0``
Performance summary workflow (targeted bed)                             ``PerformanceSummaryTargeted/v0.1.0``
Gridss                                                                  ``gridss/v2.6.2``
GATK Base Recalibration on Bam                                          ``GATKBaseRecalBQSRWorkflow/4.1.3``
GATK4 Somatic Variant Caller for Tumour Only Samples with Targeted BED  ``GATK4_SomaticVariantCallerTumorOnlyTargeted/v0.1.1``
GATK4: Haplotype Caller                                                 ``Gatk4HaplotypeCaller/4.1.3.0``
Split Multiple Alleles and Normalise Vcf                                ``SplitMultiAlleleNormaliseVcf/v0.5772``
Combine Variants                                                        ``combinevariants/0.0.8``
BGZip                                                                   ``bgzip/1.9``
BCFTools: Sort                                                          ``bcftoolssort/v1.9``
UncompressArchive                                                       ``UncompressArchive/v1.0.0``
Annotate Bam Stats to Germline Vcf Workflow                             ``AddBamStatsGermline/v0.1.0``
Tabix                                                                   ``tabix/1.2.1``
VcfLib: Vcf Length                                                      ``vcflength/v1.0.1``
VcfLib: Vcf Filter                                                      ``vcffilter/v1.0.1``
======================================================================  ======================================================



Additional configuration (inputs)
---------------------------------

======================================  ==============================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                    type                            documentation
======================================  ==============================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
sample_name                             String
fastqs                                  Array<FastqGzPair>
seqrun                                  String                          SeqRun Name (for Vcf2Tsv)
reference                               FastaWithIndexes
region_bed                              bed
region_bed_extended                     bed
region_bed_annotated                    bed
genecoverage_bed                        bed
genome_file                             TextFile
panel_name                              String
vcfcols                                 TextFile
snps_dbsnp                              CompressedIndexedVCF
snps_1000gp                             CompressedIndexedVCF
known_indels                            CompressedIndexedVCF
mills_indels                            CompressedIndexedVCF
mutalyzer_server                        String
pathos_db                               String
maxRecordsInRam                         Integer
gnomad                                  CompressedIndexedVCF
black_list                              Optional<bed>
panel_of_normals                        Optional<CompressedIndexedVCF>
fastqc_threads                          Optional<Integer>               (-t) Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine
align_and_sort_sortsam_tmpDir           Optional<String>                Undocumented option
gridss_tmpdir                           Optional<String>
haplotype_caller_pairHmmImplementation  Optional<String>                The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE
combinevariants_type                    Optional<String>                germline | somatic
combinevariants_columns                 Optional<Array<String>>         Columns to keep, seperated by space output vcf (unsorted)
filter_for_vcfs                         Optional<String>
filter_variants_1_invert                Optional<Boolean>               (-v) inverts the filter, e.g. grep -v
======================================  ==============================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================


