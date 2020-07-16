:orphan:

Annotate Bam Stats to Somatic Vcf Workflow
===============================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addBamStatsSomaticWorkflow import AddBamStatsSomatic_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addbamstatssomatic_step",
           AddBamStatsSomatic_0_1_0(
               normal_id=None,
               tumor_id=None,
               normal_bam=None,
               tumor_bam=None,
               vcf=None,
           )
       )
       wf.output("out", source=addbamstatssomatic_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for AddBamStatsSomatic:

.. code-block:: bash

   # user inputs
   janis inputs AddBamStatsSomatic > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal_bam: normal_bam.bam
       normal_id: <value>
       tumor_bam: tumor_bam.bam
       tumor_id: <value>
       vcf: vcf.vcf




5. Run AddBamStatsSomatic with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       AddBamStatsSomatic





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``AddBamStatsSomatic``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Authors: Jiaan Yu
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Embedded Tools
***************

=========================  =====================================
                           ``samtools_mpileup_subpipeline/None``
Add Bam Statistics to Vcf  ``addBamStats/0.0.7``
=========================  =====================================



Additional configuration (inputs)
---------------------------------

================  ================  ==================================
name              type              documentation
================  ================  ==================================
normal_id         String
tumor_id          String
normal_bam        IndexedBam
tumor_bam         IndexedBam
vcf               VCF
addbamstats_type  Optional<String>  must be either germline or somatic
================  ================  ==================================


