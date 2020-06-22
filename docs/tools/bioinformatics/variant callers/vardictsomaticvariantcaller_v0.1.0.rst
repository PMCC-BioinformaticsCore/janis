:orphan:

Vardict Somatic Variant Caller
============================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.vardictsomatic_variants import VardictSomaticVariantCaller

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vardictsomaticvariantcaller_step",
           VardictSomaticVariantCaller(
               normal_bam=None,
               tumor_bam=None,
               normal_name=None,
               tumor_name=None,
               intervals=None,
               header_lines=None,
               reference=None,
           )
       )
       wf.output("vardict_variants", source=vardictsomaticvariantcaller_step.vardict_variants)
       wf.output("out", source=vardictsomaticvariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vardictSomaticVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs vardictSomaticVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       header_lines: header_lines
       intervals: intervals.bed
       normal_bam: normal_bam.bam
       normal_name: <value>
       reference: reference.fasta
       tumor_bam: tumor_bam.bam
       tumor_name: <value>




5. Run vardictSomaticVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vardictSomaticVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``vardictSomaticVariantCaller``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.0
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

================  =============  ===============
name              type           documentation
================  =============  ===============
vardict_variants  CompressedVCF
out               VCF
================  =============  ===============


Embedded Tools
***************

======================  ============================
Vardict (Somatic)       ``vardict_somatic/1.6.0``
BCFTools: Annotate      ``bcftoolsAnnotate/v1.5``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
Trim IUPAC Bases        ``trimIUPAC/0.0.5``
======================  ============================



Additional configuration (inputs)
---------------------------------

============================  =================  ============================================================================
name                          type               documentation
============================  =================  ============================================================================
normal_bam                    IndexedBam
tumor_bam                     IndexedBam
normal_name                   String
tumor_name                    String
intervals                     bed
header_lines                  File
reference                     FastaWithIndexes
allele_freq_threshold         Optional<Float>
vardict_chromNamesAreNumbers  Optional<Boolean>  Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
vardict_vcfFormat             Optional<Boolean>  VCF format output
vardict_chromColumn           Optional<Integer>  The column for chromosome
vardict_regStartCol           Optional<Integer>  The column for region start, e.g. gene start
vardict_geneEndCol            Optional<Integer>  The column for region end, e.g. gene end
============================  =================  ============================================================================


