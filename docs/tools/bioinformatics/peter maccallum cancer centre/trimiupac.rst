:orphan:

Trim IUPAC Bases
============================

*0 contributors · 2 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_

Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.pmac.trimiupac.versions import TrimIUPAC_0_0_5

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "trimiupac_step",
           trimIUPAC(
               vcf=None,
           )
       )
       wf.output("out", source=trimiupac_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for trimIUPAC:

.. code-block:: bash

   # user inputs
   janis inputs trimIUPAC > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf




5. Run trimIUPAC with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       trimIUPAC





Information
------------


:ID: ``trimIUPAC``
:URL: *No URL to the documentation was provided*
:Versions: 0.0.5, 0.0.4
:Container: michaelfranklin/pmacutil:0.0.5
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  ========  ==========  ======================================
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ======================================
vcf             VCF                                    0  The VCF to remove the IUPAC bases from
outputFilename  Optional<Filename>                     2
==============  ==================  ========  ==========  ======================================
