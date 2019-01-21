
gatkhaplotypecaller
===================
*bioinformatics*

Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php# <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#/>`_

Docstring
*********
Call germline SNPs and indels via local re-assembly of haplotypes
    
    The HaplotypeCaller is capable of calling SNPs and indels simultaneously via local de-novo assembly of haplotypes 
    in an active region. In other words, whenever the program encounters a region showing signs of variation, it 
    discards the existing mapping information and completely reassembles the reads in that region. This allows the 
    HaplotypeCaller to be more accurate when calling regions that are traditionally difficult to call, for example when 
    they contain different types of variants close to each other. It also makes the HaplotypeCaller much better at 
    calling indels than position-based callers like UnifiedGenotyper.
    
    In the GVCF workflow used for scalable variant calling in DNA sequence data, HaplotypeCaller runs per-sample to 
    generate an intermediate GVCF (not to be used in final analysis), which can then be used in GenotypeGVCFs for joint 
    genotyping of multiple samples in a very efficient way. The GVCF workflow enables rapid incremental processing of 
    samples as they roll off the sequencer, as well as scaling to very large cohort sizes (e.g. the 92K exomes of ExAC).
    
    In addition, HaplotypeCaller is able to handle non-diploid organisms as well as pooled experiment data. 
    Note however that the algorithms used to calculate variant likelihoods is not well suited to extreme allele 
    frequencies (relative to ploidy) so its use is not recommended for somatic (cancer) variant discovery. 
    For that purpose, use Mutect2 instead.
    
    Finally, HaplotypeCaller is also able to correctly handle the splice junctions that make RNAseq a challenge 
    for most variant callers, on the condition that the input read data has previously been processed according 
    to our recommendations as documented (https://software.broadinstitute.org/gatk/documentation/article?id=4067).

*This page was automatically generated*
