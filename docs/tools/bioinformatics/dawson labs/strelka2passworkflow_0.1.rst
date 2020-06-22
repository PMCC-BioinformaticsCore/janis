:orphan:

Strelka 2Pass analysis
=============================================

*1 contributor · 1 version*

This is the full 2pass analysis workflow to do joint somatic variant calling with strelka2.
        The idea is similar to the RNASeq 2pass analysis, when the input of the first analysis is used to guide the second analysis.

        The workflow will
         * run manta
         * run strelka with manata output
         * run strelka with strelka and manta output
         * reannotate the filter column
         * output resuults


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.strelka2passworkflow import Strelka2PassWorkflow

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "strelka2passworkflow_step",
           Strelka2PassWorkflow(
               normalBam=None,
               tumorBams=None,
               reference=None,
           )
       )
       wf.output("snvs", source=strelka2passworkflow_step.snvs)
       wf.output("indels", source=strelka2passworkflow_step.indels)
       wf.output("svs", source=strelka2passworkflow_step.svs)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Strelka2PassWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs Strelka2PassWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normalBam: normalBam.cram
       reference: reference.fasta
       tumorBams:
       - tumorBams_0.cram
       - tumorBams_1.cram




5. Run Strelka2PassWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Strelka2PassWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``Strelka2PassWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: Sebastian Hollizeck
:Citations: 
:Created: 2019-10-11
:Updated: 2019-10-15



Outputs
-----------

======  ===========================  ===============
name    type                         documentation
======  ===========================  ===============
snvs    Array<CompressedIndexedVCF>
indels  Array<CompressedIndexedVCF>
svs     Array<CompressedIndexedVCF>
======  ===========================  ===============


Embedded Tools
***************

===============================  =================================
Strelka 2Pass analysis step1     ``Strelka2PassWorkflowStep1/0.1``
Strelka 2Pass analysis step 2    ``Strelka2PassWorkflowStep2/0.1``
Refilter Strelka2 Variant Calls  ``refilterStrelka2Calls/0.1.6``
BGZip                            ``bgzip/1.2.1``
Tabix                            ``tabix/1.2.1``
===============================  =================================



Additional configuration (inputs)
---------------------------------

===========  =======================  ===============
name         type                     documentation
===========  =======================  ===============
normalBam    CramPair
tumorBams    Array<CramPair>
reference    FastaWithIndexes
callRegions  Optional<BedTABIX>
exome        Optional<Boolean>
sampleNames  Optional<Array<String>>
===========  =======================  ===============


