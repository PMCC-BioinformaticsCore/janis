:orphan:

Filter Vardict Somatic Vcf
====================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.common.filtervardictsomaticvcf import FilterVardictSomaticVcf

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "filtervardictsomaticvcf_step",
           FilterVardictSomaticVcf(

           )
       )
       wf.output("out", source=filtervardictsomaticvcf_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FilterVardictSomaticVcf:

.. code-block:: bash

   # user inputs
   janis inputs FilterVardictSomaticVcf > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run FilterVardictSomaticVcf with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FilterVardictSomaticVcf





Information
------------


:ID: ``FilterVardictSomaticVcf``
:URL: *No URL to the documentation was provided*
:Versions: v1.9
:Container: biocontainers/bcftools:v1.9-1-deb_cv1
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

==============  ==================  ========  ==========  ===============
name            type                prefix      position  documentation
==============  ==================  ========  ==========  ===============
vcf             Optional<VCF>                          1
outputFilename  Optional<Filename>  -o                 3
==============  ==================  ========  ==========  ===============
