:orphan:

Freebayes somatic workflow
=====================================================

*0 contributors · 1 version*

This workflow uses the capabilities of freebayes to output all variants independent of the
        diploid model which then in turn allows us to create a likelihood based difference between
        the normal sample and an arbitrary amount of samples.
        This allows a joint somatic genotyping of multiple samples of the same individual.


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.dawson.workflows.freebayessomaticworkflow import FreeBayesSomaticWorkflow

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "freebayessomaticworkflow_step",
           FreeBayesSomaticWorkflow(
               bams=None,
               reference=None,
               normalSample=None,
           )
       )
       wf.output("somaticOutVcf", source=freebayessomaticworkflow_step.somaticOutVcf)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for FreeBayesSomaticWorkflow:

.. code-block:: bash

   # user inputs
   janis inputs FreeBayesSomaticWorkflow > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.cram
       - bams_1.cram
       normalSample: <value>
       reference: reference.fasta




5. Run FreeBayesSomaticWorkflow with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       FreeBayesSomaticWorkflow





Information
------------

URL: *No URL to the documentation was provided*

:ID: ``FreeBayesSomaticWorkflow``
:URL: *No URL to the documentation was provided*
:Versions: 0.1
:Authors: 
:Citations: 
:Created: 2019-10-18
:Updated: 2019-10-25



Outputs
-----------

=============  ====================  ===============
name           type                  documentation
=============  ====================  ===============
somaticOutVcf  CompressedIndexedVCF
=============  ====================  ===============


Embedded Tools
***************

====================================  ===============================
Create genomic call regions           ``CreateCallRegions/v0.1.0``
freebayes                             ``freebayes_cram/1.3.1``
Call Somatic Variants from freebayes  ``callSomaticFreeBayes/0.1.6``
VcfLib: VcfCombine                    ``vcfcombine/v1.0.1``
VcfLib: VcfStreamSort                 ``vcfstreamsort/v1.0.1``
BCFTools: Normalize                   ``bcftoolsNorm/v1.9``
VcfLib: VcfAllelicPrimitives          ``vcfallelicprimitives/v1.0.1``
VcfLib: VcfFixUp                      ``vcffixup/v1.0.1``
VcfLib: VcfUniqAlleles                ``vcfuniqalleles/v1.0.1``
VcfLib: VcfUniq                       ``vcfuniq/v1.0.1``
BGZip                                 ``bgzip/1.2.1``
Tabix                                 ``tabix/1.2.1``
====================================  ===============================



Additional configuration (inputs)
---------------------------------

================================  =======================  ================================================================================================================================================================================================================================================================================
name                              type                     documentation
================================  =======================  ================================================================================================================================================================================================================================================================================
bams                              Array<CramPair>
reference                         FastaWithIndexes
normalSample                      String
regionSize                        Optional<Integer>
sampleNames                       Optional<Array<String>>
skipCov                           Optional<Integer>
minCov                            Optional<Integer>
createCallRegions_equalize        Optional<Boolean>
callVariants_pooledDiscreteFlag   Optional<Boolean>        Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy.
callVariants_gtQuals              Optional<Boolean>        -= --genotype-qualities Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output.
callVariants_strictFlag           Optional<Boolean>        Generate strict VCF format (FORMAT/GQ will be an int)
callVariants_pooledContinousFlag  Optional<Boolean>        Output all alleles which pass input filters, regardles of genotyping outcome or model.
callVariants_reportMaxGLFlag      Optional<Boolean>        --report-genotype-likelihood-max Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods.
callVariants_noABPriorsFlag       Optional<Boolean>        -a --allele-balance-priors-off Disable use of aggregate probability of observation balance between alleles as a component of the priors.
callVariants_maxNumOfAlleles      Optional<Integer>        Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default: all)
callVariants_noPartObsFlag        Optional<Boolean>        Exclude observations which do not fully span the dynamically-determined detection window. (default, use all observations, dividing partial support across matching haplotypes when generating haplotypes.)
callVariants_useDupFlag           Optional<Boolean>        Include duplicate-marked alignments in the analysis. default: exclude duplicates marked as such in alignments
callVariants_minBaseQual          Optional<Integer>        -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q. default: 0
callVariants_minSupMQsum          Optional<Integer>        -Y --min-supporting-mapping-qsum Q Consider any allele in which and the sum of mapping qualities of supporting reads is at least Q. default: 0
callVariants_minSupQsum           Optional<Integer>        -R --min-supporting-allele-qsum Q Consider any allele in which the sum of qualities of supporting observations is at least Q. default: 0
callVariants_minAltFrac           Optional<Float>          -F --min-alternate-fraction N Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. default: 0.05
callVariants_minAltQSum           Optional<Integer>        -3 --min-alternate-qsum N Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 0
callVariants_minAltTotal          Optional<Integer>        -G --min-alternate-total N Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. default: 1
sortSomatic1_inMemoryFlag         Optional<Boolean>        load all sites and then sort in memory
normalizeSomatic1_outputType      Optional<String>         --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
normalizeSomatic1_outputFilename  Optional<Filename>       --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
allelicPrimitves_tagParsed        Optional<String>         Tag records which are split apart of a complex allele with this flag
allelicPrimitves_keepGenoFlag     Optional<Boolean>        Maintain genotype-level annotations when decomposing.  Similar caution should be used for this as for --keep-info.
sortSomatic2_inMemoryFlag         Optional<Boolean>        load all sites and then sort in memory
normalizeSomatic2_outputType      Optional<String>         --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
normalizeSomatic2_outputFilename  Optional<Filename>       --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
sortFinal_inMemoryFlag            Optional<Boolean>        load all sites and then sort in memory
================================  =======================  ================================================================================================================================================================================================================================================================================


