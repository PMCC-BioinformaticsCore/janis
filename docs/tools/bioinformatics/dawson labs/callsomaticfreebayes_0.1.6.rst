:orphan:

Call Somatic Variants from freebayes
===========================================================

*1 contributor · 1 version*

Usage: callSomaticFreeBayes.R [options]



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.callsomaticfreebayes.latest import CallSomaticFreeBayesLatest

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "callsomaticfreebayes_step",
           CallSomaticFreeBayesLatest(
               vcf=None,
           )
       )
       wf.output("out", source=callsomaticfreebayes_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for callSomaticFreeBayes:

.. code-block:: bash

   # user inputs
   janis inputs callSomaticFreeBayes > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run callSomaticFreeBayes with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       callSomaticFreeBayes





Information
------------


:ID: ``callSomaticFreeBayes``
:URL: *No URL to the documentation was provided*
:Versions: 0.1.6
:Container: shollizeck/dawsontoolkit:0.1.6.1
:Authors: Sebastian Hollizeck
:Citations: None
:Created: 2019-10-19
:Updated: 2019-10-25



Outputs
-----------

======  ======  =================
name    type    documentation
======  ======  =================
out     VCF     To determine type
======  ======  =================



Additional configuration (inputs)
---------------------------------

================  ==================  ========  ==========  ================================================================
name              type                prefix    position    documentation
================  ==================  ========  ==========  ================================================================
vcf               VCF                 -i                    input vcf
normalSampleName  Optional<String>    -n                    the normal sample name in the vcf (default: first sample in vcf)
outputFilename    Optional<Filename>  -o                    output file name (default: STDOUT)
================  ==================  ========  ==========  ================================================================
