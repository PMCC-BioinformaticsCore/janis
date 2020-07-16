:orphan:

GATK Base Recalibration on Bam
==========================================================

*0 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.gatkbasecalbam import GATKBaseRecalBQSRWorkflow_4_1_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatkbaserecalbqsrworkflow_step",
           GATKBaseRecalBQSRWorkflow_4_1_3(
               bam=None,
               reference=None,
               snps_dbsnp=None,
               snps_1000gp=None,
               known_indels=None,
               mills_indels=None,
           )
       )
       wf.output("out", source=gatkbaserecalbqsrworkflow_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for GATKBaseRecalBQSRWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs GATKBaseRecalBQSRWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       known_indels: known_indels.vcf.gz
       mills_indels: mills_indels.vcf.gz
       reference: reference.fasta
       snps_1000gp: snps_1000gp.vcf.gz
       snps_dbsnp: snps_dbsnp.vcf.gz




5. Run GATKBaseRecalBQSRWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       GATKBaseRecalBQSRWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``GATKBaseRecalBQSRWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: 4.1.2, 4.1.3
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedBam
======  ==========  ===============


Embedded Tools
***************

=============================================  =================================
GATK4: Base Recalibrator                       ``Gatk4BaseRecalibrator/4.1.3.0``
GATK4: Apply base quality score recalibration  ``Gatk4ApplyBQSR/4.1.3.0``
=============================================  =================================



Additional configuration (inputs)
---------------------------------

============  ====================  ===================================================================================================================================================
name          type                  documentation
============  ====================  ===================================================================================================================================================
bam           IndexedBam
reference     FastaWithIndexes
snps_dbsnp    CompressedIndexedVCF
snps_1000gp   CompressedIndexedVCF
known_indels  CompressedIndexedVCF
mills_indels  CompressedIndexedVCF
intervals     Optional<bed>         This optional interval supports processing by regions. If this input resolves to null, then GATK will process the whole genome per each tool's spec
============  ====================  ===================================================================================================================================================


