:orphan:

Facets: snp-pileup
====================================

*2 contributors · 3 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.facets.snp_pileup.versions import FacetsSnpPileup_0_5_14_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "facetssnppileup_step",
           FacetsSnpPileup_0_5_14_2(
               vcf_file=None,
               normal=None,
               tumour=None,
           )
       )
       wf.output("out", source=facetssnppileup_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FacetsSnpPileup:

.. code-block:: bash

   # user inputs
   janis inputs FacetsSnpPileup > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       normal: normal.bam
       tumour: tumour.bam
       vcf_file: vcf_file.vcf




5. Run FacetsSnpPileup with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FacetsSnpPileup





Information
------------


:ID: ``FacetsSnpPileup``
:URL: `https://github.com/vanallenlab/facets <https://github.com/vanallenlab/facets>`_
:Versions: 0.5.14.1, 0.5.14-2, 0.5.14
:Container: vanallenlab/facets:v0.5.14-2
:Authors: mumbler, evanwehi
:Citations: Ronglai Shen, Venkatraman E. Seshan; FACETS: allele-specific copy number and clonal heterogeneity analysis tool for high-throughput DNA sequencing, Nucleic Acids Research, Volume 44, Issue 16, 19 September 2016, Pages e131,
:DOI: https://doi.org/10.1093/nar/gkw520
:Created: 2019-12-16
:Updated: 2019-12-16



Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     File
======  ======  ===============



Additional configuration (inputs)
---------------------------------

================  ==================  ===================  ==========  =======================================================================================================
name              type                prefix                 position  documentation
================  ==================  ===================  ==========  =======================================================================================================
vcf_file          VCF                                              18
normal            IndexedBam                                       20
tumour            IndexedBam                                       21
count_orphans     Optional<Boolean>   --count-orphans               2  Do not discard anomalous read pairs
ignore_overlaps   Optional<Boolean>   --ignore-overlaps             4  Disable read-pair overlap detection.
max_depth         Optional<Integer>   --maxdepth=                   6  Sets the maximum depth. Default is 4000.
min_map_quality   Optional<Integer>   --min-map-quality=            8  Sets the minimum threshold for mapping quality. Default is 0.
min_base_quality  Optional<Integer>   --min-base-quality=          10  Sets the minimum threshold for base quality. Default is 0.
min_read_counts   Optional<String>    --min-read-counts=           12  Comma separated list of minimum read counts for a position to be output. Default is 0.
gzip              Optional<Boolean>   --gzip                       14  Compresses the output file with BGZF.
pseudo_snps       Optional<String>    --pseudo-snps=               16  Every MULTIPLE positions, if there is no SNP,insert a blank record with the total count at theposition.
output_filename   Optional<Filename>                               19
================  ==================  ===================  ==========  =======================================================================================================
