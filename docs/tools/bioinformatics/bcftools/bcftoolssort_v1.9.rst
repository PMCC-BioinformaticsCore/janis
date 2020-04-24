:orphan:

BCFTools: Sort
=============================

*0 contributors · 1 version*

About:   Sort VCF/BCF file.
Usage:   bcftools sort [OPTIONS] <FILE.vcf>


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.sort.versions import BcfToolsSort_1_9

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolssort_step",
           BcfToolsSort_1_9(
               vcf=None,
           )
       )
       wf.output("out", source=bcftoolssort_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolssort:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolssort > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run bcftoolssort with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolssort





Information
------------


:ID: ``bcftoolssort``
:URL: *No URL to the documentation was provided*
:Versions: v1.9
:Container: michaelfranklin/bcftools:1.9
:Authors: 
:Citations: None
:Created: 2019-05-09
:Updated: 2019-07-11



Outputs
-----------

======  =============  ===============
name    type           documentation
======  =============  ===============
out     CompressedVCF
======  =============  ===============



Additional configuration (inputs)
---------------------------------

==============  ==================  =============  ==========  =======================================================================================
name            type                prefix           position  documentation
==============  ==================  =============  ==========  =======================================================================================
vcf             CompressedVCF                               1  The VCF file to sort
outputFilename  Optional<Filename>  --output-file              (-o) output file name [stdout]
outputType      Optional<String>    --output-type              (-O) b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]
tempDir         Optional<String>    --temp-dir                 (-T) temporary files [/tmp/bcftools-sort.XXXXXX/]
==============  ==================  =============  ==========  =======================================================================================
