:orphan:

WGS Somatic (GATK only) [VARIANTS only]
====================================================================

``WGSSomaticGATKVariantsOnly`` · *A somatic tumor-normal variant-calling WGS pipeline using only GATK Mutect2 · 3 contributors · 1 version*

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

================  ======================  =====================================================================================================  =======================================================================================================================================================================================================================================================================================================================================================================================================================
Name              Type                    Source                                                                                                 Description
================  ======================  =====================================================================================================  =======================================================================================================================================================================================================================================================================================================================================================================================================================
gnomad            Gzipped<VCF>            * hg38: https://storage.cloud.google.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz   The genome Aggregation Database (gnomAD). This VCF must be compressed and tabix indexed. This is specific for your genome (eg: hg38 / br37) and can usually be found with your reference. For example for HG38, the Broad institute provide the following af-only-gnomad compressed and tabix indexed VCF: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=af-only
                                          * b37: https://storage.cloud.google.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf
panel_of_normals  Optional<Gzipped<VCF>>  * b37: gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf                                  VCF file of sites observed in normal.
                                          * b37-exome: gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf
reference         FastaWithIndexes        * hg38: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.fasta                     The reference genome from which to align the reads. This requires a number indexes (can be generated     with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                                                                                                                                     This pipeline expects the assembly references to be as they appear in the GCP example. For example:
                                                                                                                                                         - HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

                                                                                                                                                     - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp        Gzipped<VCF>            * hg38: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf              From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp       Gzipped<VCF>            * hg38: gs://genomics-public-data/references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz     From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``. Accessible from the HG38 genomics-public-data google cloud bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/
known_indels      Gzipped<VCF>            * hg38: gs://genomics-public-data/references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz       From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels      Gzipped<VCF>            * hg38: gs://genomics-public-data/references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz  From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
gatk_intervals    Array<bed>              None                                                                                                   List of intervals over which to split the GATK variant calling
gridss_blacklist  bed                     * hg19: https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz                  BED file containing regions to ignore. For more information, visit: https://github.com/PapenfussLab/gridss#blacklist
                                          * GRCh38: https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz
================  ======================  =====================================================================================================  =======================================================================================================================================================================================================================================================================================================================================================================================================================

4. Generate user and static input files for WGSSomaticGATKVariantsOnly:

.. code-block:: bash

   # user inputs
   janis inputs --user WGSSomaticGATKVariantsOnly > inputs.yaml

   # static inputs
   janis inputs --static WGSSomaticGATKVariantsOnly > static.yaml

**inputs.yaml**

.. code-block:: yaml

       normal_bam: NA12878-normal.bam
       normal_name: NA12878_normal
       tumor_bam: NA12878-normal.bam
       tumor_name: NA12878_tumor


**static.yaml**

.. code-block:: yaml

       gatk_intervals: BRCA1.bed
       gnomad: af-only-gnomad.hg38.vcf.gz
       gridss_blacklist: gridss_blacklist.bed
       known_indels: Homo_sapiens_assembly38.known_indels.vcf.gz
       mills_indels: Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
       reference: Homo_sapiens_assembly38.fasta
       snps_1000gp: 1000G_phase1.snps.high_confidence.hg38.vcf.gz
       snps_dbsnp: Homo_sapiens_assembly38.dbsnp138.vcf.gz


5. Run WGSSomaticGATKVariantsOnly with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       --inputs static.yaml \
       WGSSomaticGATKVariantsOnly



Outputs
-----------

===================  ============  ====================================================
name                 type          documentation
===================  ============  ====================================================
out_gridss_assembly  BAM           Assembly returned by GRIDSS
out_variants_gridss  VCF           Variants from the GRIDSS variant caller
out_variants_gatk    Gzipped<VCF>  Merged variants from the GATK caller
out_variants_split   Array<VCF>    Unmerged variants from the GATK caller (by interval)
out_variants         VCF           Final vcf
===================  ============  ====================================================

Workflow
--------

.. image:: WGSSomaticGATKVariantsOnly_1_4_0.dot.png


Information
------------


:ID: ``WGSSomaticGATKVariantsOnly``
:Versions: 1.4.0
:Authors: Michael Franklin, Richard Lupat, Jiaan Yu
:Citations: 
:Created: None
:Updated: 2020-08-18

Embedded Tools
~~~~~~~~~~~~~~~~~

==========================================  ======================================
Gridss                                      ``gridss/v2.6.2``
GATK Base Recalibration on Bam              ``GATKBaseRecalBQSRWorkflow/4.1.3``
GATK4 Somatic Variant Caller                ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs                          ``Gatk4GatherVcfs/4.1.3.0``
BGZip                                       ``bgzip/1.2.1``
BCFTools: Sort                              ``bcftoolssort/v1.9``
UncompressArchive                           ``UncompressArchive/v1.0.0``
Annotate Bam Stats to Somatic Vcf Workflow  ``AddBamStatsSomatic/v0.1.0``
==========================================  ======================================


Additional configuration (inputs)
---------------------------------

================  ======================  =======================================================================================================================================================================================================================================================================================================================================================================================================================
name              type                    documentation
================  ======================  =======================================================================================================================================================================================================================================================================================================================================================================================================================
normal_bam        IndexedBam              Indexed NORMAL bam to call somatic variants against
tumor_bam         IndexedBam              Indexed TUMOR bam to call somatic variants against
normal_name       String                  Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
tumor_name        String                  Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
gnomad            Gzipped<VCF>            The genome Aggregation Database (gnomAD). This VCF must be compressed and tabix indexed. This is specific for your genome (eg: hg38 / br37) and can usually be found with your reference. For example for HG38, the Broad institute provide the following af-only-gnomad compressed and tabix indexed VCF: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=af-only
reference         FastaWithIndexes        The reference genome from which to align the reads. This requires a number indexes (can be generated     with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

                                              This pipeline expects the assembly references to be as they appear in the GCP example. For example:
                                                  - HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

                                              - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
snps_dbsnp        Gzipped<VCF>            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
snps_1000gp       Gzipped<VCF>            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``. Accessible from the HG38 genomics-public-data google cloud bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/
known_indels      Gzipped<VCF>            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
mills_indels      Gzipped<VCF>            From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
gatk_intervals    Array<bed>              List of intervals over which to split the GATK variant calling
gridss_blacklist  bed                     BED file containing regions to ignore. For more information, visit: https://github.com/PapenfussLab/gridss#blacklist
panel_of_normals  Optional<Gzipped<VCF>>  VCF file of sites observed in normal.
================  ======================  =======================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   import "tools/gridss_v2_6_2.wdl" as G
   import "tools/GATKBaseRecalBQSRWorkflow_4_1_3.wdl" as G2
   import "tools/GATK4_SomaticVariantCaller_4_1_3_0.wdl" as G3
   import "tools/Gatk4GatherVcfs_4_1_3_0.wdl" as G4
   import "tools/bgzip_1_2_1.wdl" as B
   import "tools/bcftoolssort_v1_9.wdl" as B2
   import "tools/UncompressArchive_v1_0_0.wdl" as U
   import "tools/AddBamStatsSomatic_v0_1_0.wdl" as A

   workflow WGSSomaticGATKVariantsOnly {
     input {
       File normal_bam
       File normal_bam_bai
       File tumor_bam
       File tumor_bam_bai
       String normal_name
       String tumor_name
       File gnomad
       File gnomad_tbi
       File? panel_of_normals
       File? panel_of_normals_tbi
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       File snps_dbsnp
       File snps_dbsnp_tbi
       File snps_1000gp
       File snps_1000gp_tbi
       File known_indels
       File known_indels_tbi
       File mills_indels
       File mills_indels_tbi
       Array[File] gatk_intervals
       File gridss_blacklist
     }
     call G.gridss as vc_gridss {
       input:
         bams=[normal_bam, tumor_bam],
         bams_bai=[normal_bam_bai, tumor_bam_bai],
         reference=reference,
         reference_fai=reference_fai,
         reference_amb=reference_amb,
         reference_ann=reference_ann,
         reference_bwt=reference_bwt,
         reference_pac=reference_pac,
         reference_sa=reference_sa,
         reference_dict=reference_dict,
         blacklist=gridss_blacklist
     }
     scatter (g in gatk_intervals) {
        call G2.GATKBaseRecalBQSRWorkflow as bqsr_normal {
         input:
           bam=normal_bam,
           bam_bai=normal_bam_bai,
           intervals=g,
           reference=reference,
           reference_fai=reference_fai,
           reference_amb=reference_amb,
           reference_ann=reference_ann,
           reference_bwt=reference_bwt,
           reference_pac=reference_pac,
           reference_sa=reference_sa,
           reference_dict=reference_dict,
           snps_dbsnp=snps_dbsnp,
           snps_dbsnp_tbi=snps_dbsnp_tbi,
           snps_1000gp=snps_1000gp,
           snps_1000gp_tbi=snps_1000gp_tbi,
           known_indels=known_indels,
           known_indels_tbi=known_indels_tbi,
           mills_indels=mills_indels,
           mills_indels_tbi=mills_indels_tbi
       }
     }
     scatter (g in gatk_intervals) {
        call G2.GATKBaseRecalBQSRWorkflow as bqsr_tumor {
         input:
           bam=tumor_bam,
           bam_bai=tumor_bam_bai,
           intervals=g,
           reference=reference,
           reference_fai=reference_fai,
           reference_amb=reference_amb,
           reference_ann=reference_ann,
           reference_bwt=reference_bwt,
           reference_pac=reference_pac,
           reference_sa=reference_sa,
           reference_dict=reference_dict,
           snps_dbsnp=snps_dbsnp,
           snps_dbsnp_tbi=snps_dbsnp_tbi,
           snps_1000gp=snps_1000gp,
           snps_1000gp_tbi=snps_1000gp_tbi,
           known_indels=known_indels,
           known_indels_tbi=known_indels_tbi,
           mills_indels=mills_indels,
           mills_indels_tbi=mills_indels_tbi
       }
     }
     scatter (Q in zip(gatk_intervals, zip(transpose([bqsr_normal.out, bqsr_normal.out_bai]), transpose([bqsr_tumor.out, bqsr_tumor.out_bai])))) {
        call G3.GATK4_SomaticVariantCaller as vc_gatk {
         input:
           normal_bam=Q.right.left[0],
           normal_bam_bai=Q.right.left[1],
           tumor_bam=Q.right.right[0],
           tumor_bam_bai=Q.right.right[1],
           normal_name=normal_name,
           intervals=Q.left,
           reference=reference,
           reference_fai=reference_fai,
           reference_amb=reference_amb,
           reference_ann=reference_ann,
           reference_bwt=reference_bwt,
           reference_pac=reference_pac,
           reference_sa=reference_sa,
           reference_dict=reference_dict,
           gnomad=gnomad,
           gnomad_tbi=gnomad_tbi,
           panel_of_normals=panel_of_normals,
           panel_of_normals_tbi=panel_of_normals_tbi
       }
     }
     call G4.Gatk4GatherVcfs as vc_gatk_merge {
       input:
         vcfs=vc_gatk.out
     }
     call B.bgzip as vc_gatk_compressvcf {
       input:
         file=vc_gatk_merge.out
     }
     call B2.bcftoolssort as vc_gatk_sort_combined {
       input:
         vcf=vc_gatk_compressvcf.out
     }
     call U.UncompressArchive as vc_gatk_uncompressvcf {
       input:
         file=vc_gatk_sort_combined.out
     }
     call A.AddBamStatsSomatic as addbamstats {
       input:
         normal_id=normal_name,
         tumor_id=tumor_name,
         normal_bam=normal_bam,
         normal_bam_bai=normal_bam_bai,
         tumor_bam=tumor_bam,
         tumor_bam_bai=tumor_bam_bai,
         reference=reference,
         reference_fai=reference_fai,
         reference_amb=reference_amb,
         reference_ann=reference_ann,
         reference_bwt=reference_bwt,
         reference_pac=reference_pac,
         reference_sa=reference_sa,
         reference_dict=reference_dict,
         vcf=vc_gatk_uncompressvcf.out
     }
     output {
       File out_gridss_assembly = vc_gridss.assembly
       File out_variants_gridss = vc_gridss.out
       File out_variants_gatk = vc_gatk_sort_combined.out
       Array[File] out_variants_split = vc_gatk.out
       File out_variants = addbamstats.out
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: Workflow
   cwlVersion: v1.2
   label: WGS Somatic (GATK only) [VARIANTS only]
   doc: |
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

   requirements:
   - class: InlineJavascriptRequirement
   - class: StepInputExpressionRequirement
   - class: ScatterFeatureRequirement
   - class: SubworkflowFeatureRequirement
   - class: MultipleInputFeatureRequirement

   inputs:
   - id: normal_bam
     doc: Indexed NORMAL bam to call somatic variants against
     type: File
     secondaryFiles:
     - pattern: .bai
   - id: tumor_bam
     doc: Indexed TUMOR bam to call somatic variants against
     type: File
     secondaryFiles:
     - pattern: .bai
   - id: normal_name
     doc: |-
       Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
     type: string
   - id: tumor_name
     doc: |-
       Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
     type: string
   - id: gnomad
     doc: |-
       The genome Aggregation Database (gnomAD). This VCF must be compressed and tabix indexed. This is specific for your genome (eg: hg38 / br37) and can usually be found with your reference. For example for HG38, the Broad institute provide the following af-only-gnomad compressed and tabix indexed VCF: https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=af-only
     type: File
     secondaryFiles:
     - pattern: .tbi
   - id: panel_of_normals
     doc: VCF file of sites observed in normal.
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
   - id: reference
     doc: |2-
           The reference genome from which to align the reads. This requires a number indexes (can be generated     with the 'IndexFasta' pipeline This pipeline has been tested using the HG38 reference set.

           This pipeline expects the assembly references to be as they appear in the GCP example. For example:
               - HG38: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/

           - (".fai", ".amb", ".ann", ".bwt", ".pac", ".sa", "^.dict").
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
   - id: snps_dbsnp
     doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
     type: File
     secondaryFiles:
     - pattern: .tbi
   - id: snps_1000gp
     doc: |-
       From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``. Accessible from the HG38 genomics-public-data google cloud bucket: https://console.cloud.google.com/storage/browser/genomics-public-data/references/hg38/v0/ 
     type: File
     secondaryFiles:
     - pattern: .tbi
   - id: known_indels
     doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
     type: File
     secondaryFiles:
     - pattern: .tbi
   - id: mills_indels
     doc: From the GATK resource bundle, passed to BaseRecalibrator as ``known_sites``
     type: File
     secondaryFiles:
     - pattern: .tbi
   - id: gatk_intervals
     doc: List of intervals over which to split the GATK variant calling
     type:
       type: array
       items: File
   - id: gridss_blacklist
     doc: |-
       BED file containing regions to ignore. For more information, visit: https://github.com/PapenfussLab/gridss#blacklist
     type: File

   outputs:
   - id: out_gridss_assembly
     doc: Assembly returned by GRIDSS
     type: File
     outputSource: vc_gridss/assembly
   - id: out_variants_gridss
     doc: Variants from the GRIDSS variant caller
     type: File
     outputSource: vc_gridss/out
   - id: out_variants_gatk
     doc: Merged variants from the GATK caller
     type: File
     outputSource: vc_gatk_sort_combined/out
   - id: out_variants_split
     doc: Unmerged variants from the GATK caller (by interval)
     type:
       type: array
       items: File
     outputSource: vc_gatk/out
   - id: out_variants
     doc: Final vcf
     type: File
     outputSource: addbamstats/out

   steps:
   - id: vc_gridss
     label: Gridss
     in:
     - id: bams
       source:
       - normal_bam
       - tumor_bam
     - id: reference
       source: reference
     - id: blacklist
       source: gridss_blacklist
     run: tools/gridss_v2_6_2.cwl
     out:
     - id: out
     - id: assembly
   - id: bqsr_normal
     label: GATK Base Recalibration on Bam
     in:
     - id: bam
       source: normal_bam
     - id: intervals
       source: gatk_intervals
     - id: reference
       source: reference
     - id: snps_dbsnp
       source: snps_dbsnp
     - id: snps_1000gp
       source: snps_1000gp
     - id: known_indels
       source: known_indels
     - id: mills_indels
       source: mills_indels
     scatter:
     - intervals
     run: tools/GATKBaseRecalBQSRWorkflow_4_1_3.cwl
     out:
     - id: out
   - id: bqsr_tumor
     label: GATK Base Recalibration on Bam
     in:
     - id: bam
       source: tumor_bam
     - id: intervals
       source: gatk_intervals
     - id: reference
       source: reference
     - id: snps_dbsnp
       source: snps_dbsnp
     - id: snps_1000gp
       source: snps_1000gp
     - id: known_indels
       source: known_indels
     - id: mills_indels
       source: mills_indels
     scatter:
     - intervals
     run: tools/GATKBaseRecalBQSRWorkflow_4_1_3.cwl
     out:
     - id: out
   - id: vc_gatk
     label: GATK4 Somatic Variant Caller
     in:
     - id: normal_bam
       source: bqsr_normal/out
     - id: tumor_bam
       source: bqsr_tumor/out
     - id: normal_name
       source: normal_name
     - id: intervals
       source: gatk_intervals
     - id: reference
       source: reference
     - id: gnomad
       source: gnomad
     - id: panel_of_normals
       source: panel_of_normals
     scatter:
     - intervals
     - normal_bam
     - tumor_bam
     scatterMethod: dotproduct
     run: tools/GATK4_SomaticVariantCaller_4_1_3_0.cwl
     out:
     - id: variants
     - id: out_bam
     - id: out
   - id: vc_gatk_merge
     label: 'GATK4: Gather VCFs'
     in:
     - id: vcfs
       source: vc_gatk/out
     run: tools/Gatk4GatherVcfs_4_1_3_0.cwl
     out:
     - id: out
   - id: vc_gatk_compressvcf
     label: BGZip
     in:
     - id: file
       source: vc_gatk_merge/out
     run: tools/bgzip_1_2_1.cwl
     out:
     - id: out
   - id: vc_gatk_sort_combined
     label: 'BCFTools: Sort'
     in:
     - id: vcf
       source: vc_gatk_compressvcf/out
     run: tools/bcftoolssort_v1_9.cwl
     out:
     - id: out
   - id: vc_gatk_uncompressvcf
     label: UncompressArchive
     in:
     - id: file
       source: vc_gatk_sort_combined/out
     run: tools/UncompressArchive_v1_0_0.cwl
     out:
     - id: out
   - id: addbamstats
     label: Annotate Bam Stats to Somatic Vcf Workflow
     in:
     - id: normal_id
       source: normal_name
     - id: tumor_id
       source: tumor_name
     - id: normal_bam
       source: normal_bam
     - id: tumor_bam
       source: tumor_bam
     - id: reference
       source: reference
     - id: vcf
       source: vc_gatk_uncompressvcf/out
     run: tools/AddBamStatsSomatic_v0_1_0.cwl
     out:
     - id: out
   id: WGSSomaticGATKVariantsOnly

