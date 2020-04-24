:orphan:

Performance Validator
============================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.validation.performancevalidator import PerformanceValidator_1_2_1

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "performancevalidator_step",
           PerformanceValidator_1_2_1(
               vcf=None,
               truth=None,
               intervals=None,
           )
       )
       wf.output("summaryMetrics", source=performancevalidator_step.summaryMetrics)
   wf.output("detailMetrics", source=performancevalidator_step.detailMetrics)
   wf.output("contingencyMetrics", source=performancevalidator_step.contingencyMetrics)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for performanceValidator:

.. code-block:: bash

   # user inputs
   janis inputs performanceValidator > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       intervals:
       - intervals_0.vcf
       - intervals_1.vcf
       truth: truth.vcf
       vcf: vcf.vcf




5. Run performanceValidator with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       performanceValidator





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``performanceValidator``
:URL: *No URL to the documentation was provided*
:Versions: 1.2.1
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

==================  ======  ===============
name                type    documentation
==================  ======  ===============
summaryMetrics      File
detailMetrics       File
contingencyMetrics  File
==================  ======  ===============


Embedded Tools
***************

===========================  ====================================
BGZip                        ``bgzip/1.2.1``
Tabix                        ``tabix/1.2.1``
GATK4: Genotype Concordance  ``Gatk4GenotypeConcordance/4.1.4.0``
===========================  ====================================



Additional configuration (inputs)
---------------------------------

==========================================  =================  ========================================================================================================================================================================================================================================
name                                        type               documentation
==========================================  =================  ========================================================================================================================================================================================================================================
vcf                                         VCF
truth                                       IndexedVCF
intervals                                   Array<VCF>
genotypeConcord_treatMissingSitesAsHomeRef  Optional<Boolean>  Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth
                                                               set will be treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. This flag can only be used with a high confidence interval list.
==========================================  =================  ========================================================================================================================================================================================================================================


