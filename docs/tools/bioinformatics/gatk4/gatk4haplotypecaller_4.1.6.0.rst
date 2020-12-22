:orphan:

GATK4: Haplotype Caller
==============================================

``Gatk4HaplotypeCaller`` · *1 contributor · 7 versions*

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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.haplotypecaller.versions import Gatk4HaplotypeCaller_4_1_6

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4haplotypecaller_step",
           Gatk4HaplotypeCaller_4_1_6(
               inputRead=None,
               reference=None,
           )
       )
       wf.output("out", source=gatk4haplotypecaller_step.out)
       wf.output("bam", source=gatk4haplotypecaller_step.bam)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4HaplotypeCaller:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4HaplotypeCaller > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputRead: inputRead.bam
       reference: reference.fasta




5. Run Gatk4HaplotypeCaller with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4HaplotypeCaller





Information
------------

:ID: ``Gatk4HaplotypeCaller``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php# <https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php#>`_
:Versions: 4.1.8.1, 4.1.7.0, 4.1.6.0, 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.6.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

======  ============  ===================================================================================================
name    type          documentation
======  ============  ===================================================================================================
out     Gzipped<VCF>  A raw, unfiltered, highly sensitive callset in VCF format. File to which variants should be written
bam     IndexedBam    File to which assembled haplotypes should be written
======  ============  ===================================================================================================


Additional configuration (inputs)
---------------------------------

========================================  ========================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                                      type                      prefix                                             position  documentation
========================================  ========================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
inputRead                                 IndexedBam                --input                                                      BAM/SAM/CRAM file containing reads
reference                                 FastaWithIndexes          --reference                                               5  Reference sequence file
javaOptions                               Optional<Array<String>>
compression_level                         Optional<Integer>                                                                      Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
pairHmmImplementation                     Optional<String>          --pair-hmm-implementation                                    The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE
activityProfileOut                        Optional<String>          --activity-profile-out                                       Output the raw activity profile results in IGV format (default: null)
alleles                                   Optional<File>            --alleles                                                    (default: null) The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES
annotateWithNumDiscoveredAlleles          Optional<Boolean>         --annotate-with-num-discovered-alleles                       If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site
annotation                                Optional<Array<String>>   --annotation                                                 -A: One or more specific annotations to add to variant calls
annotationGroup                           Optional<Array<String>>   --annotation-group                                           -G	One or more groups of annotations to apply to variant calls
annotationsToExclude                      Optional<Array<String>>   --annotations-to-exclude                                     -AX	One or more specific annotations to exclude from variant calls
arguments_file                            Optional<Array<File>>     --arguments_file                                             read one or more arguments files and add them to the command line
assemblyRegionOut                         Optional<String>          --assembly-region-out                                        (default: null) Output the assembly region to this IGV formatted file. Which annotations to exclude from output in the variant calls. Note that this argument has higher priority than the -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other options.
baseQualityScoreThreshold                 Optional<Integer>         --base-quality-score-threshold                               (default: 18) Base qualities below this threshold will be reduced to the minimum (6)
cloudIndexPrefetchBuffer                  Optional<Integer>         --cloud-index-prefetch-buffer                                -CIPB (default: -1) Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.
cloudPrefetchBuffer                       Optional<Integer>         --cloud-prefetch-buffer                                      -CPB (default: 40) Size of the cloud-only prefetch buffer (in MB; 0 to disable).
contaminationFractionToFilter             Optional<Double>          --contamination-fraction-to-filter                           -contamination (default: 0.0) Fraction of contamination in sequencing data (for all samples) to aggressively remove
correctOverlappingQuality                 Optional<Boolean>         --correct-overlapping-quality                                Undocumented option
disableBamIndexCaching                    Optional<Boolean>         --disable-bam-index-caching                                  -DBIC. If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified.
founderId                                 Optional<Array<String>>   --founder-id                                                 Samples representing the population "founders"
genotypingMode                            Optional<String>          --genotyping-mode                                            (default: DISCOVERY) Specifies how to determine the alternate alleles to use for genotyping. The --genotyping-mode argument is an enumerated type (GenotypingOutputMode), which can have one of the following values: DISCOVERY (The genotyper will choose the most likely alternate allele) or GENOTYPE_GIVEN_ALLELES (Only the alleles passed by the user should be considered).
heterozygosity                            Optional<Double>          --heterozygosity                                             (default: 0.001) Heterozygosity value used to compute prior likelihoods for any locus. The expected heterozygosity value used to compute prior probability that a locus is non-reference. The default priors are for provided for humans: het = 1e-3 which means that the probability of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i) Note that heterozygosity as used here is the population genetics concept: http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics . That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. Note that this quantity has nothing to do with the likelihood of any given sample having a heterozygous genotype, which in the GATK is purely determined by the probability of the observed data P(D | AB) under the model that there may be a AB het genotype. The posterior probability of this AB genotype would use the het prior, but the GATK only uses this posterior probability in determining the prob. that a site is polymorphic. So changing the het parameters only increases the chance that a site will be called non-reference across all samples, but doesn't actually change the output genotype likelihoods at all, as these aren't posterior probabilities at all. The quantity that changes whether the GATK considers the possibility of a het genotype at all is the ploidy, which determines how many chromosomes each individual in the species carries.
heterozygosityStdev                       Optional<Double>          --heterozygosity-stdev                                       (default 0.01) Standard deviation of heterozygosity for SNP and indel calling.
indelHeterozygosity                       Optional<Double>          --indel-heterozygosity                                       (default: 1.25E-4) Heterozygosity for indel calling. This argument informs the prior probability of having an indel at a site. (See heterozygosity)
intervalMergingRule                       Optional<String>          --interval-merging-rule                                      -imr (default: ALL) Interval merging rule for abutting intervals. By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead. The --interval-merging-rule argument is an enumerated type (IntervalMergingRule), which can have one of the following values:[ALL, OVERLAPPING]
maxReadsPerAlignmentStart                 Optional<Integer>         --max-reads-per-alignment-start                              (default: 50) Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.
minBaseQualityScore                       Optional<Integer>         --min-base-quality-score                                     -mbq (default: 10) Minimum base quality required to consider a base for calling
nativePairHmmThreads                      Optional<Integer>         --native-pair-hmm-threads                                    (default: 4) How many threads should a native pairHMM implementation use
nativePairHmmUseDoublePrecision           Optional<Boolean>         --native-pair-hmm-use-double-precision                       use double precision in the native pairHmm. This is slower but matches the java implementation better
numReferenceSamplesIfNoCall               Optional<Integer>         --num-reference-samples-if-no-call                           (default: 0) Number of hom-ref genotypes to infer at sites not present in a panel. When a variant is not seen in any panel, this argument controls whether to infer (and with what effective strength) that only reference alleles were observed at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000".
outputMode                                Optional<String>          --output-mode                                                (default: EMIT_VARIANTS_ONLY) Specifies which type of calls we should output. The --output-mode argument is an enumerated type (OutputMode), which can have one of the following values: [EMIT_VARIANTS_ONLY (produces calls only at variant sites), EMIT_ALL_CONFIDENT_SITES (produces calls at variant sites and confident reference sites), EMIT_ALL_SITES (produces calls at any callable site regardless of confidence; this argument is intended only for point mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce a comprehensive set of indels in DISCOVERY mode)]
pedigree                                  Optional<File>            --pedigree                                                   -ped (default: null) Pedigree file for determining the population "founders"
populationCallset                         Optional<File>            --population-callset                                         -population (default: null) Callset to use in calculating genotype priors
sampleName                                Optional<String>          --sample-name                                                -ALIAS (default: null) Name of single sample to use from a multi-sample bam. You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This is especially useful if your samples are all in the same file but you need to run them individually through HC in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
samplePloidy                              Optional<Integer>         --sample-ploidy                                              -ploidy (default: 2) Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy). Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
sitesOnlyVcfOutput                        Optional<Boolean>         --sites-only-vcf-output                                      (default: false) If true, don't emit genotype fields when writing vcf file output.
standardMinConfidenceThresholdForCalling  Optional<Double>          --standard-min-confidence-threshold-for-calling              -stand-call-conf (default: 10.0) The minimum phred-scaled confidence threshold at which variants should be called
useNewQualCalculator                      Optional<Boolean>         --use-new-qual-calculator                                    -new-qual If provided, we will use the new AF model instead of the so-called exact model
gvcfGqBands                               Optional<Array<Integer>>  -GQB                                                         (--gvcf-gq-bands) Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)
emitRefConfidence                         Optional<String>          --emit-ref-confidence                                        (-ERC) Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)
dontUseSoftClippedBases                   Optional<Boolean>         --dont-use-soft-clipped-bases                                Do not analyze soft clipped bases in the reads
outputFilename                            Optional<Filename>        --output                                                  8  File to which variants should be written
dbsnp                                     Optional<Gzipped<VCF>>    --dbsnp                                                   7  (Also: -D) A dbSNP VCF file.
intervals                                 Optional<bed>             --intervals                                                  -L (BASE) One or more genomic intervals over which to operate
outputBamName                             Optional<Filename>        -bamout                                                   8  File to which assembled haplotypes should be written
========================================  ========================  ===============================================  ==========  =================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4HaplotypeCaller {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       String? pairHmmImplementation
       String? activityProfileOut
       File? alleles
       Boolean? annotateWithNumDiscoveredAlleles
       Array[String]? annotation
       Array[String]? annotationGroup
       Array[String]? annotationsToExclude
       Array[File]? arguments_file
       String? assemblyRegionOut
       Int? baseQualityScoreThreshold
       Int? cloudIndexPrefetchBuffer
       Int? cloudPrefetchBuffer
       Float? contaminationFractionToFilter
       Boolean? correctOverlappingQuality
       Boolean? disableBamIndexCaching
       Array[String]? founderId
       String? genotypingMode
       Float? heterozygosity
       Float? heterozygosityStdev
       Float? indelHeterozygosity
       String? intervalMergingRule
       Int? maxReadsPerAlignmentStart
       Int? minBaseQualityScore
       Int? nativePairHmmThreads
       Boolean? nativePairHmmUseDoublePrecision
       Int? numReferenceSamplesIfNoCall
       String? outputMode
       File? pedigree
       File? populationCallset
       String? sampleName
       Int? samplePloidy
       Boolean? sitesOnlyVcfOutput
       Float? standardMinConfidenceThresholdForCalling
       Boolean? useNewQualCalculator
       Array[Int]? gvcfGqBands
       String? emitRefConfidence
       Boolean? dontUseSoftClippedBases
       File inputRead
       File inputRead_bai
       File reference
       File reference_fai
       File reference_amb
       File reference_ann
       File reference_bwt
       File reference_pac
       File reference_sa
       File reference_dict
       String? outputFilename
       File? dbsnp
       File? dbsnp_tbi
       File? intervals
       String? outputBamName
     }
     command <<<
       set -e
       cp -f '~{inputRead_bai}' $(echo '~{inputRead}' | sed 's/\.[^.]*$//').bai
       gatk HaplotypeCaller \
         --java-options '-Xmx~{((select_first([runtime_memory, 8, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         ~{if defined(pairHmmImplementation) then ("--pair-hmm-implementation '" + pairHmmImplementation + "'") else ""} \
         ~{if defined(activityProfileOut) then ("--activity-profile-out '" + activityProfileOut + "'") else ""} \
         ~{if defined(alleles) then ("--alleles '" + alleles + "'") else ""} \
         ~{if (defined(annotateWithNumDiscoveredAlleles) && select_first([annotateWithNumDiscoveredAlleles])) then "--annotate-with-num-discovered-alleles" else ""} \
         ~{if (defined(annotation) && length(select_first([annotation])) > 0) then "--annotation '" + sep("' '", select_first([annotation])) + "'" else ""} \
         ~{if (defined(annotationGroup) && length(select_first([annotationGroup])) > 0) then "--annotation-group '" + sep("' '", select_first([annotationGroup])) + "'" else ""} \
         ~{if (defined(annotationsToExclude) && length(select_first([annotationsToExclude])) > 0) then "--annotations-to-exclude '" + sep("' '", select_first([annotationsToExclude])) + "'" else ""} \
         ~{if (defined(arguments_file) && length(select_first([arguments_file])) > 0) then "--arguments_file '" + sep("' '", select_first([arguments_file])) + "'" else ""} \
         ~{if defined(assemblyRegionOut) then ("--assembly-region-out '" + assemblyRegionOut + "'") else ""} \
         ~{if defined(baseQualityScoreThreshold) then ("--base-quality-score-threshold " + baseQualityScoreThreshold) else ''} \
         ~{if defined(cloudIndexPrefetchBuffer) then ("--cloud-index-prefetch-buffer " + cloudIndexPrefetchBuffer) else ''} \
         ~{if defined(cloudPrefetchBuffer) then ("--cloud-prefetch-buffer " + cloudPrefetchBuffer) else ''} \
         ~{if defined(contaminationFractionToFilter) then ("--contamination-fraction-to-filter " + contaminationFractionToFilter) else ''} \
         ~{if (defined(correctOverlappingQuality) && select_first([correctOverlappingQuality])) then "--correct-overlapping-quality" else ""} \
         ~{if (defined(disableBamIndexCaching) && select_first([disableBamIndexCaching])) then "--disable-bam-index-caching" else ""} \
         ~{if (defined(founderId) && length(select_first([founderId])) > 0) then "--founder-id '" + sep("' '", select_first([founderId])) + "'" else ""} \
         ~{if defined(genotypingMode) then ("--genotyping-mode '" + genotypingMode + "'") else ""} \
         ~{if defined(heterozygosity) then ("--heterozygosity " + heterozygosity) else ''} \
         ~{if defined(heterozygosityStdev) then ("--heterozygosity-stdev " + heterozygosityStdev) else ''} \
         ~{if defined(indelHeterozygosity) then ("--indel-heterozygosity " + indelHeterozygosity) else ''} \
         ~{if defined(intervalMergingRule) then ("--interval-merging-rule '" + intervalMergingRule + "'") else ""} \
         ~{if defined(maxReadsPerAlignmentStart) then ("--max-reads-per-alignment-start " + maxReadsPerAlignmentStart) else ''} \
         ~{if defined(minBaseQualityScore) then ("--min-base-quality-score " + minBaseQualityScore) else ''} \
         ~{if defined(nativePairHmmThreads) then ("--native-pair-hmm-threads " + nativePairHmmThreads) else ''} \
         ~{if (defined(nativePairHmmUseDoublePrecision) && select_first([nativePairHmmUseDoublePrecision])) then "--native-pair-hmm-use-double-precision" else ""} \
         ~{if defined(numReferenceSamplesIfNoCall) then ("--num-reference-samples-if-no-call " + numReferenceSamplesIfNoCall) else ''} \
         ~{if defined(outputMode) then ("--output-mode '" + outputMode + "'") else ""} \
         ~{if defined(pedigree) then ("--pedigree '" + pedigree + "'") else ""} \
         ~{if defined(populationCallset) then ("--population-callset '" + populationCallset + "'") else ""} \
         ~{if defined(sampleName) then ("--sample-name '" + sampleName + "'") else ""} \
         ~{if defined(samplePloidy) then ("--sample-ploidy " + samplePloidy) else ''} \
         ~{if (defined(sitesOnlyVcfOutput) && select_first([sitesOnlyVcfOutput])) then "--sites-only-vcf-output" else ""} \
         ~{if defined(standardMinConfidenceThresholdForCalling) then ("--standard-min-confidence-threshold-for-calling " + standardMinConfidenceThresholdForCalling) else ''} \
         ~{if (defined(useNewQualCalculator) && select_first([useNewQualCalculator])) then "--use-new-qual-calculator" else ""} \
         ~{if (defined(gvcfGqBands) && length(select_first([gvcfGqBands])) > 0) then sep(" ", prefix("-GQB ", select_first([gvcfGqBands]))) else ""} \
         ~{if defined(emitRefConfidence) then ("--emit-ref-confidence '" + emitRefConfidence + "'") else ""} \
         ~{if (defined(dontUseSoftClippedBases) && select_first([dontUseSoftClippedBases])) then "--dont-use-soft-clipped-bases" else ""} \
         --input '~{inputRead}' \
         ~{if defined(intervals) then ("--intervals '" + intervals + "'") else ""} \
         --reference '~{reference}' \
         ~{if defined(dbsnp) then ("--dbsnp '" + dbsnp + "'") else ""} \
         --output '~{select_first([outputFilename, "~{basename(inputRead, ".bam")}.vcf.gz"])}' \
         -bamout '~{select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"])}'
       if [ -f $(echo '~{select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai ]; then ln -f $(echo '~{select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"])}' | sed 's/\.[^.]*$//').bai $(echo '~{select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"])}' ).bai; fi
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.6.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "~{basename(inputRead, ".bam")}.vcf.gz"])
       File out_tbi = select_first([outputFilename, "~{basename(inputRead, ".bam")}.vcf.gz"]) + ".tbi"
       File bam = select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"])
       File bam_bai = select_first([outputBamName, "~{basename(inputRead, ".bam")}.bam"]) + ".bai"
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Haplotype Caller'
   doc: |-
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.6.0

   inputs:
   - id: javaOptions
     label: javaOptions
     type:
     - type: array
       items: string
     - 'null'
   - id: compression_level
     label: compression_level
     doc: |-
       Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
     type:
     - int
     - 'null'
   - id: pairHmmImplementation
     label: pairHmmImplementation
     doc: |-
       The PairHMM implementation to use for genotype likelihood calculations. The various implementations balance a tradeoff of accuracy and runtime. The --pair-hmm-implementation argument is an enumerated type (Implementation), which can have one of the following values: EXACT;ORIGINAL;LOGLESS_CACHING;AVX_LOGLESS_CACHING;AVX_LOGLESS_CACHING_OMP;EXPERIMENTAL_FPGA_LOGLESS_CACHING;FASTEST_AVAILABLE. Implementation:  FASTEST_AVAILABLE
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --pair-hmm-implementation
   - id: activityProfileOut
     label: activityProfileOut
     doc: 'Output the raw activity profile results in IGV format (default: null)'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --activity-profile-out
   - id: alleles
     label: alleles
     doc: |-
       (default: null) The set of alleles at which to genotype when --genotyping_mode is GENOTYPE_GIVEN_ALLELES
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --alleles
   - id: annotateWithNumDiscoveredAlleles
     label: annotateWithNumDiscoveredAlleles
     doc: |-
       If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --annotate-with-num-discovered-alleles
   - id: annotation
     label: annotation
     doc: '-A: One or more specific annotations to add to variant calls'
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --annotation
   - id: annotationGroup
     label: annotationGroup
     doc: "-G\tOne or more groups of annotations to apply to variant calls"
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --annotation-group
   - id: annotationsToExclude
     label: annotationsToExclude
     doc: "-AX\tOne or more specific annotations to exclude from variant calls"
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --annotations-to-exclude
   - id: arguments_file
     label: arguments_file
     doc: read one or more arguments files and add them to the command line
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --arguments_file
   - id: assemblyRegionOut
     label: assemblyRegionOut
     doc: |-
       (default: null) Output the assembly region to this IGV formatted file. Which annotations to exclude from output in the variant calls. Note that this argument has higher priority than the -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other options.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --assembly-region-out
   - id: baseQualityScoreThreshold
     label: baseQualityScoreThreshold
     doc: |-
       (default: 18) Base qualities below this threshold will be reduced to the minimum (6)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --base-quality-score-threshold
   - id: cloudIndexPrefetchBuffer
     label: cloudIndexPrefetchBuffer
     doc: |-
       -CIPB (default: -1) Size of the cloud-only prefetch buffer (in MB; 0 to disable). Defaults to cloudPrefetchBuffer if unset.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-index-prefetch-buffer
   - id: cloudPrefetchBuffer
     label: cloudPrefetchBuffer
     doc: '-CPB (default: 40) Size of the cloud-only prefetch buffer (in MB; 0 to disable).'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --cloud-prefetch-buffer
   - id: contaminationFractionToFilter
     label: contaminationFractionToFilter
     doc: |-
       -contamination (default: 0.0) Fraction of contamination in sequencing data (for all samples) to aggressively remove
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --contamination-fraction-to-filter
   - id: correctOverlappingQuality
     label: correctOverlappingQuality
     doc: Undocumented option
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --correct-overlapping-quality
   - id: disableBamIndexCaching
     label: disableBamIndexCaching
     doc: |-
       -DBIC. If true, don't cache bam indexes, this will reduce memory requirements but may harm performance if many intervals are specified. Caching is automatically disabled if there are no intervals specified.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --disable-bam-index-caching
   - id: founderId
     label: founderId
     doc: Samples representing the population "founders"
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --founder-id
   - id: genotypingMode
     label: genotypingMode
     doc: |-
       (default: DISCOVERY) Specifies how to determine the alternate alleles to use for genotyping. The --genotyping-mode argument is an enumerated type (GenotypingOutputMode), which can have one of the following values: DISCOVERY (The genotyper will choose the most likely alternate allele) or GENOTYPE_GIVEN_ALLELES (Only the alleles passed by the user should be considered).
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genotyping-mode
   - id: heterozygosity
     label: heterozygosity
     doc: |-
       (default: 0.001) Heterozygosity value used to compute prior likelihoods for any locus. The expected heterozygosity value used to compute prior probability that a locus is non-reference. The default priors are for provided for humans: het = 1e-3 which means that the probability of N samples being hom-ref at a site is: 1 - sum_i_2N (het / i) Note that heterozygosity as used here is the population genetics concept: http://en.wikipedia.org/wiki/Zygosity#Heterozygosity_in_population_genetics . That is, a hets value of 0.01 implies that two randomly chosen chromosomes from the population of organisms would differ from each other (one being A and the other B) at a rate of 1 in 100 bp. Note that this quantity has nothing to do with the likelihood of any given sample having a heterozygous genotype, which in the GATK is purely determined by the probability of the observed data P(D | AB) under the model that there may be a AB het genotype. The posterior probability of this AB genotype would use the het prior, but the GATK only uses this posterior probability in determining the prob. that a site is polymorphic. So changing the het parameters only increases the chance that a site will be called non-reference across all samples, but doesn't actually change the output genotype likelihoods at all, as these aren't posterior probabilities at all. The quantity that changes whether the GATK considers the possibility of a het genotype at all is the ploidy, which determines how many chromosomes each individual in the species carries.
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --heterozygosity
   - id: heterozygosityStdev
     label: heterozygosityStdev
     doc: (default 0.01) Standard deviation of heterozygosity for SNP and indel calling.
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --heterozygosity-stdev
   - id: indelHeterozygosity
     label: indelHeterozygosity
     doc: |-
       (default: 1.25E-4) Heterozygosity for indel calling. This argument informs the prior probability of having an indel at a site. (See heterozygosity)
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --indel-heterozygosity
   - id: intervalMergingRule
     label: intervalMergingRule
     doc: |-
       -imr (default: ALL) Interval merging rule for abutting intervals. By default, the program merges abutting intervals (i.e. intervals that are directly side-by-side but do not actually overlap) into a single continuous interval. However you can change this behavior if you want them to be treated as separate intervals instead. The --interval-merging-rule argument is an enumerated type (IntervalMergingRule), which can have one of the following values:[ALL, OVERLAPPING]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --interval-merging-rule
   - id: maxReadsPerAlignmentStart
     label: maxReadsPerAlignmentStart
     doc: |-
       (default: 50) Maximum number of reads to retain per alignment start position. Reads above this threshold will be downsampled. Set to 0 to disable.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-reads-per-alignment-start
   - id: minBaseQualityScore
     label: minBaseQualityScore
     doc: '-mbq (default: 10) Minimum base quality required to consider a base for calling'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --min-base-quality-score
   - id: nativePairHmmThreads
     label: nativePairHmmThreads
     doc: '(default: 4) How many threads should a native pairHMM implementation use'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --native-pair-hmm-threads
   - id: nativePairHmmUseDoublePrecision
     label: nativePairHmmUseDoublePrecision
     doc: |-
       use double precision in the native pairHmm. This is slower but matches the java implementation better
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --native-pair-hmm-use-double-precision
   - id: numReferenceSamplesIfNoCall
     label: numReferenceSamplesIfNoCall
     doc: |-
       (default: 0) Number of hom-ref genotypes to infer at sites not present in a panel. When a variant is not seen in any panel, this argument controls whether to infer (and with what effective strength) that only reference alleles were observed at that site. E.g. "If not seen in 1000Genomes, treat it as AC=0, AN=2000".
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --num-reference-samples-if-no-call
   - id: outputMode
     label: outputMode
     doc: |-
       (default: EMIT_VARIANTS_ONLY) Specifies which type of calls we should output. The --output-mode argument is an enumerated type (OutputMode), which can have one of the following values: [EMIT_VARIANTS_ONLY (produces calls only at variant sites), EMIT_ALL_CONFIDENT_SITES (produces calls at variant sites and confident reference sites), EMIT_ALL_SITES (produces calls at any callable site regardless of confidence; this argument is intended only for point mutations (SNPs) in DISCOVERY mode or generally when running in GENOTYPE_GIVEN_ALLELES mode; it will by no means produce a comprehensive set of indels in DISCOVERY mode)]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --output-mode
   - id: pedigree
     label: pedigree
     doc: '-ped (default: null) Pedigree file for determining the population "founders"'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --pedigree
   - id: populationCallset
     label: populationCallset
     doc: '-population (default: null) Callset to use in calculating genotype priors'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --population-callset
   - id: sampleName
     label: sampleName
     doc: |-
       -ALIAS (default: null) Name of single sample to use from a multi-sample bam. You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This is especially useful if your samples are all in the same file but you need to run them individually through HC in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sample-name
   - id: samplePloidy
     label: samplePloidy
     doc: |-
       -ploidy (default: 2) Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy). Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --sample-ploidy
   - id: sitesOnlyVcfOutput
     label: sitesOnlyVcfOutput
     doc: |-
       (default: false) If true, don't emit genotype fields when writing vcf file output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --sites-only-vcf-output
   - id: standardMinConfidenceThresholdForCalling
     label: standardMinConfidenceThresholdForCalling
     doc: |-
       -stand-call-conf (default: 10.0) The minimum phred-scaled confidence threshold at which variants should be called
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --standard-min-confidence-threshold-for-calling
   - id: useNewQualCalculator
     label: useNewQualCalculator
     doc: |-
       -new-qual If provided, we will use the new AF model instead of the so-called exact model
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use-new-qual-calculator
   - id: gvcfGqBands
     label: gvcfGqBands
     doc: |-
       (--gvcf-gq-bands) Exclusive upper bounds for reference confidence GQ bands (must be in [1, 100] and specified in increasing order)
     type:
     - type: array
       inputBinding:
         prefix: -GQB
       items: int
     - 'null'
     inputBinding: {}
   - id: emitRefConfidence
     label: emitRefConfidence
     doc: |-
       (-ERC) Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --emit-ref-confidence
   - id: dontUseSoftClippedBases
     label: dontUseSoftClippedBases
     doc: Do not analyze soft clipped bases in the reads
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont-use-soft-clipped-bases
   - id: inputRead
     label: inputRead
     doc: BAM/SAM/CRAM file containing reads
     type: File
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }

               return [
                       {
                           location: resolveSecondary(self.location, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     inputBinding:
       prefix: --input
   - id: reference
     label: reference
     doc: Reference sequence file
     type: File
     secondaryFiles:
     - pattern: .fai
     - pattern: .amb
     - pattern: .ann
     - pattern: .bwt
     - pattern: .pac
     - pattern: .sa
     - pattern: ^.dict
     inputBinding:
       prefix: --reference
       position: 5
   - id: outputFilename
     label: outputFilename
     doc: File to which variants should be written
     type:
     - string
     - 'null'
     default: generated.vcf.gz
     inputBinding:
       prefix: --output
       position: 8
       valueFrom: $(inputs.inputRead.basename.replace(/.bam$/, "")).vcf.gz
   - id: dbsnp
     label: dbsnp
     doc: '(Also: -D) A dbSNP VCF file.'
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --dbsnp
       position: 7
   - id: intervals
     label: intervals
     doc: -L (BASE) One or more genomic intervals over which to operate
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --intervals
   - id: outputBamName
     label: outputBamName
     doc: File to which assembled haplotypes should be written
     type:
     - string
     - 'null'
     default: generated.bam
     inputBinding:
       prefix: -bamout
       position: 8
       valueFrom: $(inputs.inputRead.basename.replace(/.bam$/, "")).bam

   outputs:
   - id: out
     label: out
     doc: |-
       A raw, unfiltered, highly sensitive callset in VCF format. File to which variants should be written
     type: File
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $(inputs.inputRead.basename.replace(/.bam$/, "")).vcf.gz
       loadContents: false
   - id: bam
     label: bam
     doc: File to which assembled haplotypes should be written
     type: File
     secondaryFiles:
     - |-
       ${

               function resolveSecondary(base, secPattern) {
                 if (secPattern[0] == "^") {
                   var spl = base.split(".");
                   var endIndex = spl.length > 1 ? spl.length - 1 : 1;
                   return resolveSecondary(spl.slice(undefined, endIndex).join("."), secPattern.slice(1));
                 }
                 return base + secPattern
               }
               return [
                       {
                           path: resolveSecondary(self.path, "^.bai"),
                           basename: resolveSecondary(self.basename, ".bai"),
                           class: "File",
                       }
               ];

       }
     outputBinding:
       glob: $(inputs.inputRead.basename.replace(/.bam$/, "")).bam
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - HaplotypeCaller
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 8, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4HaplotypeCaller


