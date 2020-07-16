:orphan:

Molpath Germline Workflow
===================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.molpathGermlineWorkflow import MolpathGermline_1_0_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "molpathgermlineworkflow_step",
           MolpathGermline_1_0_0(
               sample_name=None,
               fastqs=None,
               reference=None,
               region_bed=None,
               region_bed_extended=None,
               region_bed_annotated=None,
               genecoverage_bed=None,
               genome_file=None,
               snps_dbsnp=None,
               snps_1000gp=None,
               known_indels=None,
               mills_indels=None,
           )
       )
       wf.output("fastq_qc", source=molpathgermlineworkflow_step.fastq_qc)
       wf.output("markdups_bam", source=molpathgermlineworkflow_step.markdups_bam)
       wf.output("doc", source=molpathgermlineworkflow_step.doc)
       wf.output("summary", source=molpathgermlineworkflow_step.summary)
       wf.output("gene_summary", source=molpathgermlineworkflow_step.gene_summary)
       wf.output("region_summary", source=molpathgermlineworkflow_step.region_summary)
       wf.output("gridss_vcf", source=molpathgermlineworkflow_step.gridss_vcf)
       wf.output("gridss_bam", source=molpathgermlineworkflow_step.gridss_bam)
       wf.output("hap_vcf", source=molpathgermlineworkflow_step.hap_vcf)
       wf.output("hap_bam", source=molpathgermlineworkflow_step.hap_bam)
       wf.output("normalise_vcf", source=molpathgermlineworkflow_step.normalise_vcf)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for MolpathGermlineWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs MolpathGermlineWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       fastqs:
       - - fastqs_0.fastq.gz
         - fastqs_1.fastq.gz
       - - fastqs_0.fastq.gz
         - fastqs_1.fastq.gz
       genecoverage_bed: genecoverage_bed.bed
       genome_file: genome_file.txt
       known_indels: known_indels.vcf.gz
       mills_indels: mills_indels.vcf.gz
       reference: reference.fasta
       region_bed: region_bed.bed
       region_bed_annotated: region_bed_annotated.bed
       region_bed_extended: region_bed_extended.bed
       sample_name: <value>
       snps_1000gp: snps_1000gp.vcf.gz
       snps_dbsnp: snps_dbsnp.vcf.gz




5. Run MolpathGermlineWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       MolpathGermlineWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``MolpathGermlineWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: v1.0.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

==============  ====================  ===============
name            type                  documentation
==============  ====================  ===============
fastq_qc        Array<Array<Zip>>
markdups_bam    IndexedBam
doc             TextFile
summary         csv
gene_summary    TextFile
region_summary  TextFile
gridss_vcf      VCF
gridss_bam      BAM
hap_vcf         CompressedIndexedVCF
hap_bam         IndexedBam
normalise_vcf   VCF
==============  ====================  ===============


Embedded Tools
***************

===========================================  ========================================
FastQC                                       ``fastqc/v0.11.5``
Parse FastQC Adaptors                        ``ParseFastqcAdaptors/v0.1.0``
Align and sort reads                         ``BwaAligner/1.0.0``
Merge and Mark Duplicates                    ``mergeAndMarkBams/4.1.3``
Annotate GATK3 DepthOfCoverage Workflow      ``AnnotateDepthOfCoverage/v0.1.0``
Performance summary workflow (targeted bed)  ``PerformanceSummaryTargeted/v0.1.0``
Gridss                                       ``gridss/v2.6.2``
GATK Base Recalibration on Bam               ``GATKBaseRecalBQSRWorkflow/4.1.3``
GATK4: Haplotype Caller                      ``Gatk4HaplotypeCaller/4.1.3.0``
Split Multiple Alleles and Normalise Vcf     ``SplitMultiAlleleNormaliseVcf/v0.5772``
Annotate Bam Stats to Germline Vcf Workflow  ``AddBamStatsGermline/v0.1.0``
===========================================  ========================================



Additional configuration (inputs)
---------------------------------

======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                    type                  documentation
======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================
sample_name                             String
fastqs                                  Array<FastqGzPair>
reference                               FastaWithIndexes
region_bed                              bed
region_bed_extended                     bed
region_bed_annotated                    bed
genecoverage_bed                        bed
genome_file                             TextFile
snps_dbsnp                              CompressedIndexedVCF
snps_1000gp                             CompressedIndexedVCF
known_indels                            CompressedIndexedVCF
mills_indels                            CompressedIndexedVCF
black_list                              Optional<bed>
fastqc_threads                          Optional<Integer>     (-t) Specifies the number of files which can be processed simultaneously. Each thread will be allocated 250MB of memory so you shouldn't run more threads than your available memory will cope with, and not more than 6 threads on a 32 bit machine
align_and_sort_sortsam_tmpDir           Optional<String>      Undocumented option
gridss_tmpdir                           Optional<String>
haplotype_caller_pairHmmImplementation  Optional<String>      The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE
======================================  ====================  =============================================================================================================================================================================================================================================================================================================================================================================================================================================


