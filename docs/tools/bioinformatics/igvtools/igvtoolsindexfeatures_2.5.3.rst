:orphan:

IGVTools: Index Features
================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.igvtools.index.versions import IgvIndexFeature_2_5_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "igvtoolsindexfeatures_step",
           IgvIndexFeature_2_5_3(
               inp=None,
           )
       )
       wf.output("out", source=igvtoolsindexfeatures_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for IgvToolsIndexFeatures:

.. code-block:: bash

   # user inputs
   janis inputs IgvToolsIndexFeatures > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inp: inp.vcf




5. Run IgvToolsIndexFeatures with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       IgvToolsIndexFeatures





Information
------------


:ID: ``IgvToolsIndexFeatures``
:URL: *No URL to the documentation was provided*
:Versions: 2.5.3
:Container: quay.io/biocontainers/igvtools:2.5.3--0
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ==========  ===============
name    type        documentation
======  ==========  ===============
out     IndexedVCF
======  ==========  ===============



Additional configuration (inputs)
---------------------------------

======  ======  ========  ==========  ===============
name    type    prefix      position  documentation
======  ======  ========  ==========  ===============
inp     VCF                        1
======  ======  ========  ==========  ===============
