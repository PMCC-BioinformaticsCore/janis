
gatk4genotypeconcordance
========================
*bioinformatics*

Documentation
-------------

URL
******
`https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php>`_

Docstring
*********
GenotypeConcordance (Picard)
    Calculates the concordance between genotype data of one samples in each of two VCFs - one being 
    considered the truth (or reference) the other being the call. The concordance is broken into 
    separate results sections for SNPs and indels. Statistics are reported in three different files.
    
    Summary
    Calculates the concordance between genotype data of one samples in each of two VCFs - one being 
    considered the truth (or reference) the other being the call. The concordance is broken into 
    separate results sections for SNPs and indels. Summary and detailed statistics are reported.
    
    Details
    This tool evaluates the concordance between genotype calls for a sample in different callsets
    where one is being considered as the "truth" (aka standard, or reference) and the other as the 
    "call" that is being evaluated for accuracy. The Comparison can be restricted to a confidence 
    interval which is typically used in order to enable proper assessment of False Positives and 
    the False-Positive Rate (FPR).
     
    Output Metrics:
    Output metrics consists of GenotypeConcordanceContingencyMetrics, GenotypeConcordanceSummaryMetrics, 
    and GenotypeConcordanceDetailMetrics. For each set of metrics, the data is broken into separate 
    sections for SNPs and INDELs. Note that only SNP and INDEL variants are considered, MNP, Symbolic, 
    and Mixed classes of variants are not included.
    
    GenotypeConcordanceContingencyMetrics enumerate the constituents of each contingent in a callset 
    including true-positive (TP), true-negative (TN), false-positive (FP), and false-negative (FN) calls.
    GenotypeConcordanceDetailMetrics include the numbers of SNPs and INDELs for each contingent genotype 
    as well as the number of validated genotypes.
    
    GenotypeConcordanceSummaryMetrics provide specific details for the variant caller performance 
    on a callset including values for sensitivity, specificity, and positive predictive values.
    
    
    Useful definitions applicable to alleles and genotypes:
        - Truthset - A callset (typically in VCF format) containing variant calls and genotypes that have been 
            cross-validated with multiple technologies e.g. Genome In A Bottle Consortium (GIAB) (https://sites.stanford.edu/abms/giab)
        - TP - True-positives are variant sites that match against the truth-set
        - FP - False-positives are reference sites miscalled as variant
        - FN - False-negatives are variant sites miscalled as reference
        - TN - True-negatives are correctly called as reference
        - Validated genotypes - are TP sites where the exact genotype (HET or HOM-VAR) appears in the truth-set
    
    VCF Output:
        - The concordance state will be stored in the CONC_ST tag in the INFO field
        - The truth sample name will be "truth" and call sample name will be "call"

*This page was automatically generated*
