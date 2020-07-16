:orphan:

Vardict Germline Variant Caller
==============================================================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.variantcallers.vardictgermline_variants import VardictGermlineVariantCaller

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vardictgermlinevariantcaller_step",
           VardictGermlineVariantCaller(
               bam=None,
               intervals=None,
               sample_name=None,
               header_lines=None,
               reference=None,
           )
       )
       wf.output("variants", source=vardictgermlinevariantcaller_step.variants)
       wf.output("out", source=vardictgermlinevariantcaller_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vardictGermlineVariantCaller:

.. code-block:: bash

   # user inputs
   janis inputs vardictGermlineVariantCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       header_lines: header_lines
       intervals: intervals.bed
       reference: reference.fasta
       sample_name: <value>




5. Run vardictGermlineVariantCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vardictGermlineVariantCaller





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``vardictGermlineVariantCaller``
:URL: *No URL to the documentation was provided*
:Versions: v0.1.1
:Authors: 
:Citations: 
:Created: None
:Updated: None



Outputs
-----------

========  ====================  ===============
name      type                  documentation
========  ====================  ===============
variants  CompressedIndexedVCF
out       VCF
========  ====================  ===============


Embedded Tools
***************

======================  ============================
VarDict (Germline)      ``vardict_germline/1.6.0``
BCFTools: Annotate      ``bcftoolsAnnotate/v1.5``
BGZip                   ``bgzip/1.2.1``
Tabix                   ``tabix/1.2.1``
Split Multiple Alleles  ``SplitMultiAllele/v0.5772``
Trim IUPAC Bases        ``trimIUPAC/0.0.5``
VcfTools                ``VcfTools/0.1.16``
======================  ============================



Additional configuration (inputs)
---------------------------------

=============================  =================  =================================================================================================================================================================================================================================================================
name                           type               documentation
=============================  =================  =================================================================================================================================================================================================================================================================
bam                            IndexedBam
intervals                      bed
sample_name                    String
header_lines                   File
reference                      FastaWithIndexes
allele_freq_threshold          Optional<Float>
vardict_chromNamesAreNumbers   Optional<Boolean>  Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
vardict_vcfFormat              Optional<Boolean>  VCF format output
vardict_chromColumn            Optional<Integer>  The column for chromosome
vardict_regStartCol            Optional<Integer>  The column for region start, e.g. gene start
vardict_geneEndCol             Optional<Integer>  The column for region end, e.g. gene end
compressvcf_stdout             Optional<Boolean>  c: Write to standard output, keep original files unchanged.
filterpass_removeFileteredAll  Optional<Boolean>  Removes all sites with a FILTER flag other than PASS.
filterpass_recode              Optional<Boolean>
filterpass_recodeINFOAll       Optional<Boolean>  These options can be used with the above recode options to define an INFO key name to keep in the output  file.  This  option can be used multiple times to keep more of the INFO fields. The second option is used to keep all INFO values in the original file.
=============================  =================  =================================================================================================================================================================================================================================================================


