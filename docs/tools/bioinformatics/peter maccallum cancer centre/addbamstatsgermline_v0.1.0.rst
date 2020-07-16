:orphan:

Annotate Bam Stats to Germline Vcf Workflow
=================================================================

*1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.addBamStatsGermlineWorkflow import AddBamStatsGermline_0_1_0

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "addbamstatsgermline_step",
           AddBamStatsGermline_0_1_0(
               bam=None,
               vcf=None,
           )
       )
       wf.output("out", source=addbamstatsgermline_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for AddBamStatsGermline:

.. code-block:: bash

   # user inputs
   janis inputs AddBamStatsGermline > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       vcf: vcf.vcf




5. Run AddBamStatsGermline with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       AddBamStatsGermline





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``AddBamStatsGermline``
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

=========================  =========================
SamTools: Mpileup          ``SamToolsMpileup/1.9.0``
Add Bam Statistics to Vcf  ``addBamStats/0.0.7``
=========================  =========================



Additional configuration (inputs)
---------------------------------

============================  =================  ========================================================
name                          type               documentation
============================  =================  ========================================================
bam                           IndexedBam
vcf                           VCF
samtoolsmpileup_countOrphans  Optional<Boolean>  do not discard anomalous read pairs
samtoolsmpileup_noBAQ         Optional<Boolean>  disable BAQ (per-Base Alignment Quality)
samtoolsmpileup_minBQ         Optional<Integer>  Minimum base quality for a base to be considered [13]
samtoolsmpileup_maxDepth      Optional<Integer>  max per-file depth; avoids excessive memory usage [8000]
addbamstats_type              Optional<String>   must be either germline or somatic
============================  =================  ========================================================


