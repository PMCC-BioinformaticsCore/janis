:orphan:

VcfTools: VcfConcat
=======================================

*1 contributor · 1 version*

Concatenates VCF files (for example split by chromosome). Note that the input and output VCFs will have the same number of columns, the script does not merge VCFs by position (see also vcf-merge).

In the basic mode it does not do anything fancy except for a sanity check that all files have the same columns. When run with the -s option, it will perform a partial merge sort, looking at limited number of open files simultaneously.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcftools.vcfconcat.versions import VcfToolsVcfConcat_0_1_16

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcftoolsvcfconcat_step",
           VcfToolsVcfConcat_0_1_16(
               vcfTabix=None,
           )
       )
       wf.output("out", source=vcftoolsvcfconcat_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for VcfToolsVcfConcat:

.. code-block:: bash

   # user inputs
   janis inputs VcfToolsVcfConcat > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcfTabix:
       - vcfTabix_0.vcf.gz
       - vcfTabix_1.vcf.gz




5. Run VcfToolsVcfConcat with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       VcfToolsVcfConcat





Information
------------


:ID: ``VcfToolsVcfConcat``
:URL: `http://vcftools.sourceforge.net/perl_module.html#vcf-concat <http://vcftools.sourceforge.net/perl_module.html#vcf-concat>`_
:Versions: 0.1.16
:Container: biocontainers/vcftools:v0.1.16-1-deb_cv1
:Authors: Jiaan Yu
:Citations: None
:Created: 2020-05-21
:Updated: 2020-05-21



Outputs
-----------

======  ===========  ===============
name    type         documentation
======  ===========  ===============
out     stdout<VCF>
======  ===========  ===============



Additional configuration (inputs)
---------------------------------

============  ===========================  ============  ==========  =============================================================================
name          type                         prefix          position  documentation
============  ===========================  ============  ==========  =============================================================================
vcfTabix      Array<CompressedIndexedVCF>                        10
checkColumns  Optional<Boolean>            -c                        Do not concatenate, only check if the columns agree.
padMissing    Optional<Boolean>            -p                        Write '.' in place of missing columns. Useful for joining chrY with the rest.
mergeSort     Optional<Integer>            --merge-sort              Allow small overlaps in N consecutive files.
============  ===========================  ============  ==========  =============================================================================
