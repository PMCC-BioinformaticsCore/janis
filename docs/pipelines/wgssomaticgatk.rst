:orphan:

WGS Somatic (GATK only)
========================================

A somatic tumor-normal variant-calling WGS pipeline using only GATK Mutect2 · 1 contributor · 1 version

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Quickstart
-----------

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate an inputs file for WGSSomaticGATK:

.. code-block:: bash
   
   janis inputs WGSSomaticGATK > inputs.yaml

**inputs.yaml**

.. code-block:: yaml

       gatk_intervals:
       - gatk_intervals_0.bed
       - gatk_intervals_1.bed
       known_indels: known_indels.vcf.gz
       mills_indels: mills_indels.vcf.gz
       normal_inputs:
       - - normal_inputs_0.fastq.gz
         - normal_inputs_1.fastq.gz
       - - normal_inputs_0.fastq.gz
         - normal_inputs_1.fastq.gz
       normal_name: <value>
       reference: reference.fasta
       snps_1000gp: snps_1000gp.vcf.gz
       snps_dbsnp: snps_dbsnp.vcf.gz
       tumor_inputs:
       - - tumor_inputs_0.fastq.gz
         - tumor_inputs_1.fastq.gz
       - - tumor_inputs_0.fastq.gz
         - tumor_inputs_1.fastq.gz
       tumor_name: <value>


5. Run the WGSSomaticGATK pipeline with:

.. code-block:: bash

   janis run [...workflow options] --inputs inputs.yaml WGSSomaticGATK



Outputs
-----------

=============  =================  ===============
name           type               documentation
=============  =================  ===============
normal_bam     IndexedBam
tumor_bam      IndexedBam
normal_report  Array<Array<Zip>>
tumor_report   Array<Array<Zip>>
variants_gatk  CompressedVCF
=============  =================  ===============


Information
------------

:ID: ``WGSSomaticGATK``
:Python: ``janis_pipelines.wgs_somatic_gatk.wgssomaticgatk import WGSSomaticGATK``
:Versions: 1.2.0
:Authors: Michael Franklin
:Citations: 
:Created: None
:Updated: 2019-10-16

Embedded Tools
~~~~~~~~~~~~~~~~~

============================  ======================================================================================================================================
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10852e358>>``
                              ``somatic_subpipeline/<bound method WorkflowBuilder.version of <janis_core.workflow.workflow.WorkflowBuilder object at 0x10852d0f0>>``
GATK4 Somatic Variant Caller  ``GATK4_SomaticVariantCaller/4.1.3.0``
GATK4: Gather VCFs            ``Gatk4GatherVcfs/4.1.3.0``
BCFTools: Sort                ``bcftoolssort/v1.9``
============================  ======================================================================================================================================


Additional configuration (inputs)
---------------------------------

=================  ====================  =======================================================================================================================================================================================================================================================================
name               type                  documentation
=================  ====================  =======================================================================================================================================================================================================================================================================
normal_inputs      Array<FastqGzPair>    An array of NORMAL FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
tumor_inputs       Array<FastqGzPair>    An array of TUMOR FastqGz pairs. These are aligned separately and merged to create higher depth coverages from multiple sets of reads
normal_name        String                Sample name for the NORMAL sample from which to generate the readGroupHeaderLine for BwaMem
tumor_name         String                Sample name for the TUMOR sample from which to generate the readGroupHeaderLine for BwaMem
gatk_intervals     Array<bed>            List of intervals over which to split the GATK variant calling
reference          FastaWithIndexes      The reference genome from which to align the reads. This requires a number indexes (can be generated with the 'IndexFasta' pipeline. This pipeline has been tested with the hg38 reference genome.
snps_dbsnp         CompressedIndexedVCF  From the GATK resource bundle
snps_1000gp        CompressedIndexedVCF  From the GATK resource bundle
known_indels       CompressedIndexedVCF  From the GATK resource bundle
mills_indels       CompressedIndexedVCF  From the GATK resource bundle
cutadapt_adapters  Optional<File>        Specifies a file which contains a list of sequences to determine valid overrepresented sequences from the FastQC report to trim with Cuatadapt. The file must contain sets of named adapters in the form name[tab]sequence. Lines prefixed with a hash will be ignored.
=================  ====================  =======================================================================================================================================================================================================================================================================
