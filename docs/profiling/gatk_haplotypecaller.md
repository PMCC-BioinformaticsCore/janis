# GATK4: Haplotype Caller

Haplotype Caller is a GATK tool for calling SNPs and indels. Its 
[documentation page](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php)
lists how the tool works. This page specifically talks about the GATK4 version, 
though it's broadly applicable to lower versions as well.

Optimising Haplotype Caller isn't intuitive, and simply giving the tool 
more resources does not have predictable results. For this reason, a lot of the GATK tools disable
multi-threading in favour of the software: [Spark](#).

Specifically, without Sparks, the [HaplotypeCaller will not utilise any multithreaded operation](https://gatkforums.broadinstitute.org/gatk/discussion/comment/54982/#Comment_54982).

, GATK (together with Intel) have produced a separate:

- **Beta**: [HaplotypeCallerSpark](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_HaplotypeCallerSpark.php)