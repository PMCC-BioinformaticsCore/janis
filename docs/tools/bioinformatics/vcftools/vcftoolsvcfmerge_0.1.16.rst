:orphan:

VcfTools: VcfMerge
=====================================

*1 contributor · 1 version*

Merges two or more VCF files into one so that, for example, if two source files had one column each, on output will be printed a file with two columns. See also vcf-concat for concatenating VCFs split by chromosome.

vcf-merge A.vcf.gz B.vcf.gz C.vcf.gz | bgzip -c > out.vcf.gz

Note that this script is not intended for concatenating VCF files. For this, use vcf-concat instead.
Note: A fast htslib C version of this tool is now available (see bcftools merge).


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vcftools.vcfmerge.versions import VcfToolsVcfMerge_0_1_16

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vcftoolsvcfmerge_step",
           VcfToolsVcfMerge_0_1_16(
               vcfTabix=None,
           )
       )
       wf.output("out", source=vcftoolsvcfmerge_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for VcfToolsVcfMerge:

.. code-block:: bash

   # user inputs
   janis inputs VcfToolsVcfMerge > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcfTabix:
       - vcfTabix_0.vcf.gz
       - vcfTabix_1.vcf.gz




5. Run VcfToolsVcfMerge with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       VcfToolsVcfMerge





Information
------------


:ID: ``VcfToolsVcfMerge``
:URL: `http://vcftools.sourceforge.net/perl_module.html#vcf-merge <http://vcftools.sourceforge.net/perl_module.html#vcf-merge>`_
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

================  ===========================  ===================  ==========  ===================================================================================================================================================================
name              type                         prefix                 position  documentation
================  ===========================  ===================  ==========  ===================================================================================================================================================================
vcfTabix          Array<CompressedIndexedVCF>                               10
collapse          Optional<String>             -c                               treat as identical sites with differing alleles [any] <snps|indels|both|any|none>
removeDuplicates  Optional<Boolean>            --remove-duplicates              If there should be two consecutive rows with the same chr:pos, print only the first one.
vcfHeader         Optional<File>               --vcf-header                     Use the provided VCF header
regionsList       Optional<Array<String>>      --regions                        Do only the given regions (comma-separated list).
regionsFile       Optional<File>               --regions                        Do only the given regions (one region per line in a file).
refForMissing     Optional<String>             --ref-for-missing                Use the REF allele instead of the default missing genotype. Because it is not obvious what ploidy should be used, a user-defined string is used instead (e.g. 0/0).
silent            Optional<Boolean>            --silent                         Try to be a bit more silent, no warnings about duplicate lines.
trimALTs          Optional<Boolean>            --trim-ALTs                      If set, redundant ALTs will be removed
================  ===========================  ===================  ==========  ===================================================================================================================================================================
