:orphan:

WGS Somatic (Multi callers)
====================================================

*A somatic tumor-normal variant-calling WGS pipeline using GATK, VarDict and Strelka2 · 3 contributors · 1 version*

This is a genomics pipeline to align sequencing data (Fastq pairs) into BAMs:

- Takes raw sequence data in the FASTQ format;
- align to the reference genome using BWA MEM;
- Marks duplicates using Picard;
- Call the appropriate somatic variant callers (GATK / Strelka / VarDict);
- Outputs the final variants in the VCF format.

**Resources**

This pipeline has been tested using the HG38 reference set, available on Google Cloud Storage through:

- https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

This pipeline expects the assembly references to be as they appear in that storage     (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
The known sites (snps_dbsnp, snps_1000gp, known_indels, mills_indels) should be gzipped and tabix indexed.


Quickstart
-----------

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.

====================  ====================  =========================================================================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================
Name                  Type                  Example                                                                                                                                                                                    Description
====================  ====================  =========================================================================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================
cutadapt_adapters     Optional<File>        https://github.com/csf-ngs/fastqc/blob/master/Contaminants/contaminant_list.txt                                                                                                            Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
gatk_intervals        Array<bed>            BRCA1.bed                                                                                                                                                                                  List of intervals over which to split the GATK variant calling
gridss_blacklist      bed                   https://github.com/PapenfussLab/gridss#blacklist                                                                                                                                           BED file containing regions to ignore.
vardict_intervals     Array<bed>            BRCA1.bed                                                                                                                                                                                  List of intervals over which to split the VarDict variant calling
strelka_intervals     BedTABIX              BRCA1.bed.gz                                                                                                                                                                               An interval for which to restrict the analysis to.
vardict_header_lines  File                  https://gist.githubusercontent.com/illusional/5b75a0506f7327aca7d355f8ad5008f8/raw/e181c0569771e6a557d01a8a1f70c71e3598a269/headerLines.txt                                                As with chromosomal sequences it is highly recommended (but not required) that the header include tags describing the contigs referred to in the VCF file. This furthermore allows these contigs to come from different files. The format is identical to that of a reference sequence, but with an additional URL tag to indicate where that sequence can be found. For example:

                                                                                                                                                                                                                                       .. code-block:

                                                                                                                                                                                                                                          ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

                                                                                                                                                                                                                                       Source: (1.2.5 Alternative allele field format) https://samtools.github.io/hts-specs/VCFv4.1.pdf (edited)
reference             FastaWithIndexes      HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta                                                                                                           This pipeline expects the assembly references to be as they appear in the GCP example:

                                                                                                                                                                                                                                       - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp            CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            (WARNING: The file available from the genomics-public-data resource on Google Cloud Storage is NOT compressed and indexed. This will need to be completed prior to starting the pipeline.

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
snps_1000gp           CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
known_indels          CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
mills_indels          CompressedIndexedVCF  HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/                                                                                            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``

                                            File: gs://genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
====================  ====================  =========================================================================================================================================================================================  =================================================================================================================================================================================================================================================================================================================================================================================

4. Generate user and static input files for WGSSomaticMultiCallers:

.. code-block:: bash

   # user inputs
   janis inputs --user WGSSomaticMultiCallers > inputs.yaml

   # static inputs
   janis inputs --static WGSSomaticMultiCallers > static.yaml

**inputs.yaml**

.. code-block:: yaml

       normal_inputs:
       - - normal_R1.fastq.gz
         - normal_R2.fastq.gz
       - - normal_R1-TOPUP.fastq.gz
         - normal_R2-TOPUP.fastq.gz
       normal_name: <value>
       tumor_inputs:
       - - tumor_R1.fastq.gz
         - tumor_R2.fastq.gz
       - - tumor_R1-TOPUP.fastq.gz
         - tumor_R2-TOPUP.fastq.gz
       tumor_name: <value>


**static.yaml**

.. code-block:: yaml

       gatk_intervals:
       - gatk_intervals_0.bed
       - gatk_intervals_1.bed
       gridss_blacklist: gridss_blacklist.bed
       known_indels: Homo_sapiens_assembly38.known_indels.vcf.gz
       mills_indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
       reference: Homo_sapiens_assembly38.fasta
       snps_1000gp: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
       snps_dbsnp: Homo_sapiens_assembly38.dbsnp138.vcf.gz
       strelka_intervals: strelka_intervals.bed.gz
       vardict_header_lines: vardict_header_lines
       vardict_intervals:
       - vardict_intervals_0.bed
       - vardict_intervals_1.bed


5. Run WGSSomaticMultiCallers with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --inputs static.yaml \
       WGSSomaticMultiCallers



Outputs
-----------

================  =================  ================================================
name              type               documentation
================  =================  ================================================
normal_report     Array<Array<Zip>>  A zip file of the NORMAL FastQC quality reports.
tumor_report      Array<Array<Zip>>  A zip file of the TUMOR FastQC quality reports.
normal_bam        IndexedBam         Aligned and indexed NORMAL bam
tumor_bam         IndexedBam         Aligned and indexed TUMOR bam
gridss_assembly   VCF                Assembly returned by GRIDSS
variants_gatk     VCF                Merged variants from the GATK caller
variants_strelka  VCF                Variants from the Strelka variant caller
variants_vardict  VCF                Merged variants from the VarDict caller
variants_gridss   VCF                Variants from the GRIDSS variant caller
variants          VCF                Combined variants from all 3 callers
================  =================  ================================================


Information
------------

:ID: ``WGSSomaticMultiCallers``
:Versions: 1.2.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: 2018-12-24
:Updated: 2020-03-05

Embedded Tools
~~~~~~~~~~~~~~~~~

==============================  ======================================================================================================================================
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x1093397b8>>``
                                ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x109352ba8>>``
GATK4 Somatic Variant Caller    ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs              ``Gatk4GatherVcfs/4.1.3.0``
Strelka Somatic Variant Caller  ``strelkaSomaticVariantCaller/v0.1.0``
Gridss                          ``gridss/v2.5.1-dev``
Vardict Somatic Variant Caller  ``vardictSomaticVariantCaller/v0.1.0``
Combine Variants                ``combinevariants/0.0.4``
BCFTools: Sort                  ``bcftoolssort/v1.9``
==============================  ======================================================================================================================================


Additional configuration (inputs)
---------------------------------

========================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
name                      type                     documentation
========================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
normal_inputs             Array<FastqGzPair>       An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
tumor_inputs              Array<FastqGzPair>       An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
normal_name               String                   Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
tumor_name                String                   Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
gatk_intervals            Array<bed>               List of intervals over which to split the GATK variant calling
gridss_blacklist          bed                      BED file containing regions to ignore.
vardict_intervals         Array<bed>               List of intervals over which to split the VarDict variant calling
strelka_intervals         BedTABIX                 An interval for which to restrict the analysis to.
vardict_header_lines      File                     As with chromosomal sequences it is highly recommended (but not required) that the header include tags describing the contigs referred to in the VCF file. This furthermore allows these contigs to come from different files. The format is identical to that of a reference sequence, but with an additional URL tag to indicate where that sequence can be found. For example:

                                                   .. code-block:

                                                      ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>

                                                   Source: (1.2.5 Alternative allele field format) https://samtools.github.io/hts-specs/VCFv4.1.pdf (edited)
reference                 FastaWithIndexes         The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                                   This pipeline expects the assembly references to be as they appear in the GCP example:

                                                   - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp                CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp               CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
known_indels              CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels              CompressedIndexedVCF     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
cutadapt_adapters         Optional<File>           Specifies a containment list for cutadapt, which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form: ``name[tab]sequence``. Lines prefixed with a hash will be ignored.
allele_freq_threshold     Optional<Float>          The threshold for VarDict's allele frequency, default: 0.05 or 5%
combine_variants_type     Optional<String>         germline | somatic
combine_variants_columns  Optional<Array<String>>  Columns to keep, seperated by space output vcf (unsorted)
========================  =======================  =================================================================================================================================================================================================================================================================================================================================================================================
