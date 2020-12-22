:orphan:

freebayes
=========

``freebayes`` · *2 contributors · 2 versions*

usage: freebayes [OPTION] ... [BAM FILE] ...
Bayesian haplotype-based polymorphism discovery.
Version:1.2.0



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.freebayes.versions import FreeBayes_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "freebayes_step",
           FreeBayes_1_2(
               bams=None,
               reference=None,
               theta=None,
               ploidy=None,
               refQual=None,
               maxNumOfAlleles=None,
               haplotypeLength=None,
               minRepSize=None,
               minRepEntropy=None,
               useDupFlag=None,
               minMappingQual=None,
               minBaseQual=None,
               minSupQsum=None,
               minSupMQsum=None,
               minSupBQthres=None,
               maxMisMatchFrac=None,
               minAltFrac=None,
               minAltCount=None,
               minAltQSum=None,
               minAltTotal=None,
               minCov=None,
               probContamin=None,
               genotypingMaxIter=None,
               genotypingMaxBDepth=None,
               postIntegrationLim=None,
               readDepFact=None,
           )
       )
       wf.output("out", source=freebayes_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for freebayes:

.. code-block:: bash

   # user inputs
   janis inputs freebayes > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bams:
       - bams_0.bam
       - bams_1.bam
       reference: reference.fasta




5. Run freebayes with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       freebayes





Information
------------

:ID: ``freebayes``
:URL: `https://github.com/ekg/freebayes <https://github.com/ekg/freebayes>`_
:Versions: 1.3.1, 1.2
:Container: papaemmelab/docker-freebayes:v0.1.5
:Authors: Sebastian Hollizeck, Michael Franklin
:Citations: Garrison E, Marth G. Haplotype-based variant detection from short-read sequencing. arXiv preprint arXiv:1207.3907 [q-bio.GN] 2012
:Created: 2019-10-08
:Updated: 2019-10-19


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================
name                   type                prefix                            position    documentation
=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================
bams                   Array<IndexedBam>   -b                                            Add FILE to the set of BAM files to be analyzed.
reference              FastaFai            -f                                            Use FILE as the reference sequence for analysis. An index file (FILE.fai) will be created if none exists. If neither --targets nor --region are specified, FreeBayes will analyze every position in this reference.
theta                  Float               -T                                            The expected mutation rate or pairwise nucleotide diversity among the population under analysis. This serves as the single parameter to the Ewens Sampling Formula prior model default: 0.001
ploidy                 Integer             -p                                            Sets the default ploidy for the analysis to N. default: 2
refQual                String              --reference-quality                           --reference-quality MQ,BQ  Assign mapping quality of MQ to the reference allele at each site and base quality of BQ. default: 100,60
maxNumOfAlleles        Integer             -n                                            Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default: all)
haplotypeLength        Integer             --haplotype-length                            Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping. (default: 3)
minRepSize             Integer             --min-repeat-size                             When assembling observations across repeats, require the total repeat length at least this many bp. (default: 5)
minRepEntropy          Integer             --min-repeat-entropy                          To detect interrupted repeats, build across sequence until it has  entropy > N bits per bp. Set to 0 to turn off. (default: 1)
useDupFlag             Boolean             -4                                            Include duplicate-marked alignments in the analysis. default: exclude duplicates marked as such in alignments
minMappingQual         Integer             -m                                            Exclude alignments from analysis if they have a mapping quality less than Q. default: 1
minBaseQual            Integer             -q                                            -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q. default: 0
minSupQsum             Integer             -R                                            -R --min-supporting-allele-qsum Q Consider any allele in which the sum of qualities of supporting observations is at least Q. default: 0
minSupMQsum            Integer             -Y                                            -Y --min-supporting-mapping-qsum Q Consider any allele in which and the sum of mapping qualities of supporting reads is at least Q. default: 0
minSupBQthres          Integer             -Q                                            -Q --mismatch-base-quality-threshold Q Count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= Q. default: 10
maxMisMatchFrac        Float               -z                                            -z --read-max-mismatch-fraction N Exclude reads with more than N [0,1] fraction of mismatches where each mismatch has base quality >= mismatch-base-quality-threshold default: 1.0
minAltFrac             Float               -F                                            -F --min-alternate-fraction N Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. default: 0.05
minAltCount            Integer             -C                                            -C --min-alternate-count N Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 2
minAltQSum             Integer             -3                                            -3 --min-alternate-qsum N Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 0
minAltTotal            Integer             -G                                            -G --min-alternate-total N Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. default: 1
minCov                 Integer             --min-coverage                                --min-coverage N Require at least this coverage to process a site. default: 0
probContamin           Float               --prob-contamination                          --prob-contamination F An estimate of contamination to use for all samples. default: 10e-9
genotypingMaxIter      Integer             -B                                            -B --genotyping-max-iterations N Iterate no more than N times during genotyping step. default: 1000.
genotypingMaxBDepth    Integer             --genotyping-max-banddepth                    --genotyping-max-banddepth N Integrate no deeper than the Nth best genotype by likelihood when genotyping. default: 6.
postIntegrationLim     String              -W                                            -W --posterior-integration-limits N,M Integrate all genotype combinations in our posterior space which include no more than N samples with their Mth best data likelihood. default: 1,3.
readDepFact            Float               -D                                            -D --read-dependence-factor N Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations. default: 0.9
bamList                Optional<TextFile>  -L                                            A file containing a list of BAM files to be analyzed.
targetsFile            Optional<bed>       -t                                            Limit analysis to targets listed in the BED-format FILE.
region                 Optional<String>    -r                                            <chrom>:<start_position>-<end_position> Limit analysis to the specified region, 0-base coordinates, end_position not included (same as BED format). Either '-' or '..' maybe used as a separator.
samplesFile            Optional<TextFile>  -s                                            FILE  Limit analysis to samples listed (one per line) in the FILE. By default FreeBayes will analyze all samples in its input BAM files.
popFile                Optional<TextFile>  --populations                                 FILE Each line of FILE should list a sample and a population which it is part of. The population-based bayesian inference model will then be partitioned on the basis of the populations.
cnvFile                Optional<TextFile>  -A                                            FILE Read a copy number map from the BED file FILE, which has either a sample-level ploidy: sample name, copy number or a region-specific format: reference sequence, start, end, sample name, copy number ... for each region in each sample which does not have the default copy number as set by --ploidy.
outputFilename         Optional<Filename>  -v                                            FILE Output VCF-format results to FILE. (default: stdout)
gvcfFlag               Optional<Boolean>   --gvcf                                        Write gVCF output, which indicates coverage in uncalled regions.
gvcfChunkSize          Optional<Integer>   --gvcf-chunk                                  When writing gVCF output emit a record for every NUM bases.
candidateVcf           Optional<File>      -@                                            Use variants reported in VCF file as input to the algorithm. Variants in this file will included in the output even if there is not enough support in the data to pass input filters.
restrictSitesFlag      Optional<Boolean>   -l                                            Only provide variant calls and genotype likelihoods for sites and alleles which are provided in the VCF input, and provide output in the VCF for all input alleles, not just those which have support in the data.
candidateHaploVcf      Optional<File>      --haplotype-basis-alleles                     When specified, only variant alleles provided in this input VCF will be used for the construction of complex or haplotype alleles.
reportHapAllelesFlag   Optional<Boolean>   --report-all-haplotype-alleles                At sites where genotypes are made over haplotype alleles, provide information about all alleles in output, not only those which are called.
monomorphicFlag        Optional<Boolean>   --report-monomorphic                          Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Loci which do not have any potential alternates have '.' for ALT.
polyMoprhProbFlag      Optional<Float>     -P                                            Report sites if the probability that there is a polymorphism at the site is greater than N. default: 0.0. Note that post-filtering is generally recommended over the use of this parameter.
strictFlag             Optional<Boolean>   --strict-vcf                                  Generate strict VCF format (FORMAT/GQ will be an int)
pooledDiscreteFlag     Optional<Boolean>   -J                                            Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy.
pooledContinousFlag    Optional<Boolean>   -K                                            Output all alleles which pass input filters, regardles of genotyping outcome or model.
addRefFlag             Optional<Boolean>   -Z                                            This flag includes the reference allele in the analysis as if it is another sample from the same population.
ignoreSNPsFlag         Optional<Boolean>   -I                                            Ignore SNP alleles.
ignoreINDELsFlag       Optional<Boolean>   -i                                            Ignore insertion and deletion alleles.
ignoreMNPsFlag         Optional<Boolean>   -X                                            Ignore multi-nuceotide polymorphisms, MNPs.
ignoreComplexVarsFlag  Optional<Boolean>   -u                                            Ignore complex events (composites of other classes).
maxNumOfComplexVars    Optional<Integer>   -E
noPartObsFlag          Optional<Boolean>   --no-partial-observations                     Exclude observations which do not fully span the dynamically-determined detection window. (default, use all observations, dividing partial support across matching haplotypes when generating haplotypes.)
noNormaliseFlag        Optional<Boolean>   -O                                            Turn off left-alignment of indels, which is enabled by default.
readMisMatchLim        Optional<Integer>   -U                                            -U --read-mismatch-limit N Exclude reads with more than N mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. default: ~unbounded
readSNPLim             Optional<Integer>   -$                                            -$ --read-snp-limit N Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. default: ~unbounded
readINDELLim           Optional<Integer>   -e                                            -e --read-indel-limit N Exclude reads with more than N separate gaps. default: ~unbounded
standardFilterFlag     Optional<Boolean>   -0                                            -0 --standard-filters Use stringent input base and mapping quality filters Equivalent to -m 30 -q 20 -R 0 -S 0
maxCov                 Optional<Integer>   --max-coverage                                --max-coverage N Do not process sites with greater than this coverage. default: no limit
noPopPriorsFlag        Optional<Boolean>   -k                                            -k --no-population-priors Equivalent to --pooled-discrete --hwe-priors-off and removal of Ewens Sampling Formula component of priors.
noHWEPriorsFlag        Optional<Boolean>   -w                                            -w --hwe-priors-off Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by observation frequency.
noBinOBSPriorsFlag     Optional<Boolean>   -V                                            -V --binomial-obs-priors-off Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and read position (5'-3') probability.
noABPriorsFlag         Optional<Boolean>   -a                                            -a --allele-balance-priors-off Disable use of aggregate probability of observation balance between alleles as a component of the priors.
obsBiasFile            Optional<TextFile>  --observation-bias                            --observation-bias FILE Read length-dependent allele observation biases from FILE. The format is [length] [alignment efficiency relative to reference] where the efficiency is 1 if there is no relative observation bias.
baseQualCap            Optional<Integer>   --base-quality-cap                            --base-quality-cap Q Limit estimated observation quality by capping base quality at Q.
legGLScalc             Optional<Boolean>   --legacy-gls                                  --legacy-gls Use legacy (polybayes equivalent) genotype likelihood calculations
contaminEst            Optional<TextFile>  --contamination-estimates                     --contamination-estimates FILE A file containing per-sample estimates of contamination, such as those generated by VerifyBamID. The format should be: sample p(read=R|genotype=AR) p(read=A|genotype=AA) Sample '*' can be used to set default contamination estimates.
repoprtMaxGLFlag       Optional<Boolean>   --report-genotype-likelihood-max              --report-genotype-likelihood-max Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods.
excludeUnObsGT         Optional<Boolean>   -N                                            -N --exclude-unobserved-genotypes Skip sample genotypings for which the sample has no supporting reads.
gtVarThres             Optional<Integer>   -S                                            -S --genotype-variant-threshold N Limit posterior integration to samples where the second-best genotype likelihood is no more than log(N) from the highest genotype likelihood for the sample. default: ~unbounded
useMQFlag              Optional<Boolean>   -j                                            -j --use-mapping-quality Use mapping quality of alleles when calculating data likelihoods.
harmIndelQualFlag      Optional<Boolean>   -H                                            -H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel. By default use a minimum BQ in flanking sequence.
gtQuals                Optional<Boolean>   -=                                            -= --genotype-qualities Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output.
=====================  ==================  ================================  ==========  =============================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task freebayes {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[File] bams
       Array[File] bams_bai
       File? bamList
       File reference
       File reference_fai
       File? targetsFile
       String? region
       File? samplesFile
       File? popFile
       File? cnvFile
       String? outputFilename
       Boolean? gvcfFlag
       Int? gvcfChunkSize
       File? candidateVcf
       Boolean? restrictSitesFlag
       File? candidateHaploVcf
       Boolean? reportHapAllelesFlag
       Boolean? monomorphicFlag
       Float? polyMoprhProbFlag
       Boolean? strictFlag
       Float? theta
       Int? ploidy
       Boolean? pooledDiscreteFlag
       Boolean? pooledContinousFlag
       Boolean? addRefFlag
       String? refQual
       Boolean? ignoreSNPsFlag
       Boolean? ignoreINDELsFlag
       Boolean? ignoreMNPsFlag
       Boolean? ignoreComplexVarsFlag
       Int? maxNumOfAlleles
       Int? maxNumOfComplexVars
       Int? haplotypeLength
       Int? minRepSize
       Int? minRepEntropy
       Boolean? noPartObsFlag
       Boolean? noNormaliseFlag
       Boolean? useDupFlag
       Int? minMappingQual
       Int? minBaseQual
       Int? minSupQsum
       Int? minSupMQsum
       Int? minSupBQthres
       Int? readMisMatchLim
       Float? maxMisMatchFrac
       Int? readSNPLim
       Int? readINDELLim
       Boolean? standardFilterFlag
       Float? minAltFrac
       Int? minAltCount
       Int? minAltQSum
       Int? minAltTotal
       Int? minCov
       Int? maxCov
       Boolean? noPopPriorsFlag
       Boolean? noHWEPriorsFlag
       Boolean? noBinOBSPriorsFlag
       Boolean? noABPriorsFlag
       File? obsBiasFile
       Int? baseQualCap
       Float? probContamin
       Boolean? legGLScalc
       File? contaminEst
       Boolean? repoprtMaxGLFlag
       Int? genotypingMaxIter
       Int? genotypingMaxBDepth
       String? postIntegrationLim
       Boolean? excludeUnObsGT
       Int? gtVarThres
       Boolean? useMQFlag
       Boolean? harmIndelQualFlag
       Float? readDepFact
       Boolean? gtQuals
     }
     command <<<
       set -e
       freebayes \
         ~{if length(bams) > 0 then "-b '" + sep("' -b '", bams) + "'" else ""} \
         ~{if defined(bamList) then ("-L '" + bamList + "'") else ""} \
         -f '~{reference}' \
         ~{if defined(targetsFile) then ("-t '" + targetsFile + "'") else ""} \
         ~{if defined(region) then ("-r '" + region + "'") else ""} \
         ~{if defined(samplesFile) then ("-s '" + samplesFile + "'") else ""} \
         ~{if defined(popFile) then ("--populations '" + popFile + "'") else ""} \
         ~{if defined(cnvFile) then ("-A '" + cnvFile + "'") else ""} \
         -v '~{select_first([outputFilename, "generated.vcf"])}' \
         ~{if select_first([gvcfFlag, false]) then "--gvcf" else ""} \
         ~{if defined(gvcfChunkSize) then ("--gvcf-chunk " + gvcfChunkSize) else ''} \
         ~{if defined(candidateVcf) then ("-@ '" + candidateVcf + "'") else ""} \
         ~{if (defined(restrictSitesFlag) && select_first([restrictSitesFlag])) then "-l" else ""} \
         ~{if defined(candidateHaploVcf) then ("--haplotype-basis-alleles '" + candidateHaploVcf + "'") else ""} \
         ~{if (defined(reportHapAllelesFlag) && select_first([reportHapAllelesFlag])) then "--report-all-haplotype-alleles" else ""} \
         ~{if (defined(monomorphicFlag) && select_first([monomorphicFlag])) then "--report-monomorphic" else ""} \
         ~{if defined(select_first([polyMoprhProbFlag, 0.0])) then ("-P " + select_first([polyMoprhProbFlag, 0.0])) else ''} \
         ~{if (defined(strictFlag) && select_first([strictFlag])) then "--strict-vcf" else ""} \
         -T ~{select_first([theta, 0.001])} \
         -p ~{select_first([ploidy, 2])} \
         ~{if (defined(pooledDiscreteFlag) && select_first([pooledDiscreteFlag])) then "-J" else ""} \
         ~{if (defined(pooledContinousFlag) && select_first([pooledContinousFlag])) then "-K" else ""} \
         ~{if (defined(addRefFlag) && select_first([addRefFlag])) then "-Z" else ""} \
         --reference-quality '~{select_first([refQual, "100,60"])}' \
         ~{if (defined(ignoreSNPsFlag) && select_first([ignoreSNPsFlag])) then "-I" else ""} \
         ~{if (defined(ignoreINDELsFlag) && select_first([ignoreINDELsFlag])) then "-i" else ""} \
         ~{if (defined(ignoreMNPsFlag) && select_first([ignoreMNPsFlag])) then "-X" else ""} \
         ~{if (defined(ignoreComplexVarsFlag) && select_first([ignoreComplexVarsFlag])) then "-u" else ""} \
         -n ~{select_first([maxNumOfAlleles, 0])} \
         ~{if defined(maxNumOfComplexVars) then ("-E " + maxNumOfComplexVars) else ''} \
         --haplotype-length ~{select_first([haplotypeLength, 3])} \
         --min-repeat-size ~{select_first([minRepSize, 5])} \
         --min-repeat-entropy ~{select_first([minRepEntropy, 1])} \
         ~{if (defined(noPartObsFlag) && select_first([noPartObsFlag])) then "--no-partial-observations" else ""} \
         ~{if (defined(noNormaliseFlag) && select_first([noNormaliseFlag])) then "-O" else ""} \
         ~{if select_first([useDupFlag, false]) then "-4" else ""} \
         -m ~{select_first([minMappingQual, 1])} \
         -q ~{select_first([minBaseQual, 0])} \
         -R ~{select_first([minSupQsum, 0])} \
         -Y ~{select_first([minSupMQsum, 0])} \
         -Q ~{select_first([minSupBQthres, 10])} \
         ~{if defined(readMisMatchLim) then ("-U " + readMisMatchLim) else ''} \
         -z ~{select_first([maxMisMatchFrac, 1.0])} \
         ~{if defined(readSNPLim) then ("-$ " + readSNPLim) else ''} \
         ~{if defined(readINDELLim) then ("-e " + readINDELLim) else ''} \
         ~{if (defined(standardFilterFlag) && select_first([standardFilterFlag])) then "-0" else ""} \
         -F ~{select_first([minAltFrac, 0.05])} \
         -C ~{select_first([minAltCount, 2])} \
         -3 ~{select_first([minAltQSum, 0])} \
         -G ~{select_first([minAltTotal, 1])} \
         --min-coverage ~{select_first([minCov, 0])} \
         ~{if defined(maxCov) then ("--max-coverage " + maxCov) else ''} \
         ~{if (defined(noPopPriorsFlag) && select_first([noPopPriorsFlag])) then "-k" else ""} \
         ~{if (defined(noHWEPriorsFlag) && select_first([noHWEPriorsFlag])) then "-w" else ""} \
         ~{if (defined(noBinOBSPriorsFlag) && select_first([noBinOBSPriorsFlag])) then "-V" else ""} \
         ~{if (defined(noABPriorsFlag) && select_first([noABPriorsFlag])) then "-a" else ""} \
         ~{if defined(obsBiasFile) then ("--observation-bias '" + obsBiasFile + "'") else ""} \
         ~{if defined(baseQualCap) then ("--base-quality-cap " + baseQualCap) else ''} \
         --prob-contamination ~{select_first([probContamin, 1e-09])} \
         ~{if (defined(legGLScalc) && select_first([legGLScalc])) then "--legacy-gls" else ""} \
         ~{if defined(contaminEst) then ("--contamination-estimates '" + contaminEst + "'") else ""} \
         ~{if (defined(repoprtMaxGLFlag) && select_first([repoprtMaxGLFlag])) then "--report-genotype-likelihood-max" else ""} \
         -B ~{select_first([genotypingMaxIter, 1000])} \
         --genotyping-max-banddepth ~{select_first([genotypingMaxBDepth, 6])} \
         -W '~{select_first([postIntegrationLim, "1,3"])}' \
         ~{if (defined(excludeUnObsGT) && select_first([excludeUnObsGT])) then "-N" else ""} \
         ~{if defined(gtVarThres) then ("-S " + gtVarThres) else ''} \
         ~{if (defined(useMQFlag) && select_first([useMQFlag])) then "-j" else ""} \
         ~{if (defined(harmIndelQualFlag) && select_first([harmIndelQualFlag])) then "-H" else ""} \
         -D ~{select_first([readDepFact, 0.9])} \
         ~{if (defined(gtQuals) && select_first([gtQuals])) then "-=" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "papaemmelab/docker-freebayes:v0.1.5"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 16, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: freebayes
   doc: |
     usage: freebayes [OPTION] ... [BAM FILE] ...
     Bayesian haplotype-based polymorphism discovery.
     Version:1.2.0

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: papaemmelab/docker-freebayes:v0.1.5

   inputs:
   - id: bams
     label: bams
     doc: Add FILE to the set of BAM files to be analyzed.
     type:
       type: array
       inputBinding:
         prefix: -b
       items: File
     inputBinding: {}
   - id: bamList
     label: bamList
     doc: A file containing a list of BAM files to be analyzed.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -L
   - id: reference
     label: reference
     doc: |2-
        Use FILE as the reference sequence for analysis. An index file (FILE.fai) will be created if none exists. If neither --targets nor --region are specified, FreeBayes will analyze every position in this reference.
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: -f
   - id: targetsFile
     label: targetsFile
     doc: ' Limit analysis to targets listed in the BED-format FILE.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -t
   - id: region
     label: region
     doc: |-
       <chrom>:<start_position>-<end_position> Limit analysis to the specified region, 0-base coordinates, end_position not included (same as BED format). Either '-' or '..' maybe used as a separator.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -r
   - id: samplesFile
     label: samplesFile
     doc: |-
       FILE  Limit analysis to samples listed (one per line) in the FILE. By default FreeBayes will analyze all samples in its input BAM files.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -s
   - id: popFile
     label: popFile
     doc: |-
       FILE Each line of FILE should list a sample and a population which it is part of. The population-based bayesian inference model will then be partitioned on the basis of the populations.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --populations
   - id: cnvFile
     label: cnvFile
     doc: |-
       FILE Read a copy number map from the BED file FILE, which has either a sample-level ploidy: sample name, copy number or a region-specific format: reference sequence, start, end, sample name, copy number ... for each region in each sample which does not have the default copy number as set by --ploidy.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -A
   - id: outputFilename
     label: outputFilename
     doc: 'FILE Output VCF-format results to FILE. (default: stdout)'
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: -v
   - id: gvcfFlag
     label: gvcfFlag
     doc: Write gVCF output, which indicates coverage in uncalled regions.
     type: boolean
     default: false
     inputBinding:
       prefix: --gvcf
   - id: gvcfChunkSize
     label: gvcfChunkSize
     doc: ' When writing gVCF output emit a record for every NUM bases.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --gvcf-chunk
   - id: candidateVcf
     label: candidateVcf
     doc: |2-
        Use variants reported in VCF file as input to the algorithm. Variants in this file will included in the output even if there is not enough support in the data to pass input filters.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -@
   - id: restrictSitesFlag
     label: restrictSitesFlag
     doc: |-
       Only provide variant calls and genotype likelihoods for sites and alleles which are provided in the VCF input, and provide output in the VCF for all input alleles, not just those which have support in the data.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -l
   - id: candidateHaploVcf
     label: candidateHaploVcf
     doc: |-
       When specified, only variant alleles provided in this input VCF will be used for the construction of complex or haplotype alleles.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --haplotype-basis-alleles
   - id: reportHapAllelesFlag
     label: reportHapAllelesFlag
     doc: |-
       At sites where genotypes are made over haplotype alleles, provide information about all alleles in output, not only those which are called.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --report-all-haplotype-alleles
   - id: monomorphicFlag
     label: monomorphicFlag
     doc: |2-
        Report even loci which appear to be monomorphic, and report all considered alleles, even those which are not in called genotypes. Loci which do not have any potential alternates have '.' for ALT.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --report-monomorphic
   - id: polyMoprhProbFlag
     label: polyMoprhProbFlag
     doc: |-
       Report sites if the probability that there is a polymorphism at the site is greater than N. default: 0.0. Note that post-filtering is generally recommended over the use of this parameter.
     type: float
     default: 0.0
     inputBinding:
       prefix: -P
   - id: strictFlag
     label: strictFlag
     doc: Generate strict VCF format (FORMAT/GQ will be an int)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --strict-vcf
   - id: theta
     label: theta
     doc: |-
       The expected mutation rate or pairwise nucleotide diversity among the population under analysis. This serves as the single parameter to the Ewens Sampling Formula prior model default: 0.001
     type: float
     default: 0.001
     inputBinding:
       prefix: -T
   - id: ploidy
     label: ploidy
     doc: 'Sets the default ploidy for the analysis to N. default: 2'
     type: int
     default: 2
     inputBinding:
       prefix: -p
   - id: pooledDiscreteFlag
     label: pooledDiscreteFlag
     doc: |-
       Assume that samples result from pooled sequencing. Model pooled samples using discrete genotypes across pools. When using this flag, set --ploidy to the number of alleles in each sample or use the --cnv-map to define per-sample ploidy.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -J
   - id: pooledContinousFlag
     label: pooledContinousFlag
     doc: |-
       Output all alleles which pass input filters, regardles of genotyping outcome or model.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -K
   - id: addRefFlag
     label: addRefFlag
     doc: |-
       This flag includes the reference allele in the analysis as if it is another sample from the same population.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -Z
   - id: refQual
     label: refQual
     doc: |-
       --reference-quality MQ,BQ  Assign mapping quality of MQ to the reference allele at each site and base quality of BQ. default: 100,60
     type: string
     default: 100,60
     inputBinding:
       prefix: --reference-quality
   - id: ignoreSNPsFlag
     label: ignoreSNPsFlag
     doc: Ignore SNP alleles.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -I
   - id: ignoreINDELsFlag
     label: ignoreINDELsFlag
     doc: Ignore insertion and deletion alleles.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -i
   - id: ignoreMNPsFlag
     label: ignoreMNPsFlag
     doc: Ignore multi-nuceotide polymorphisms, MNPs.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -X
   - id: ignoreComplexVarsFlag
     label: ignoreComplexVarsFlag
     doc: Ignore complex events (composites of other classes).
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -u
   - id: maxNumOfAlleles
     label: maxNumOfAlleles
     doc: |-
       Evaluate only the best N SNP alleles, ranked by sum of supporting quality scores. (Set to 0 to use all; default: all)
     type: int
     default: 0
     inputBinding:
       prefix: -n
   - id: maxNumOfComplexVars
     label: maxNumOfComplexVars
     doc: ''
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -E
   - id: haplotypeLength
     label: haplotypeLength
     doc: |-
       Allow haplotype calls with contiguous embedded matches of up to this length. Set N=-1 to disable clumping. (default: 3)
     type: int
     default: 3
     inputBinding:
       prefix: --haplotype-length
   - id: minRepSize
     label: minRepSize
     doc: |-
       When assembling observations across repeats, require the total repeat length at least this many bp. (default: 5)
     type: int
     default: 5
     inputBinding:
       prefix: --min-repeat-size
   - id: minRepEntropy
     label: minRepEntropy
     doc: |-
       To detect interrupted repeats, build across sequence until it has  entropy > N bits per bp. Set to 0 to turn off. (default: 1)
     type: int
     default: 1
     inputBinding:
       prefix: --min-repeat-entropy
   - id: noPartObsFlag
     label: noPartObsFlag
     doc: |-
       Exclude observations which do not fully span the dynamically-determined detection window. (default, use all observations, dividing partial support across matching haplotypes when generating haplotypes.)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-partial-observations
   - id: noNormaliseFlag
     label: noNormaliseFlag
     doc: Turn off left-alignment of indels, which is enabled by default.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -O
   - id: useDupFlag
     label: useDupFlag
     doc: |-
       Include duplicate-marked alignments in the analysis. default: exclude duplicates marked as such in alignments
     type: boolean
     default: false
     inputBinding:
       prefix: '-4'
   - id: minMappingQual
     label: minMappingQual
     doc: |2-
        Exclude alignments from analysis if they have a mapping quality less than Q. default: 1
     type: int
     default: 1
     inputBinding:
       prefix: -m
   - id: minBaseQual
     label: minBaseQual
     doc: |2-
        -q --min-base-quality Q Exclude alleles from analysis if their supporting base quality is less than Q. default: 0
     type: int
     default: 0
     inputBinding:
       prefix: -q
   - id: minSupQsum
     label: minSupQsum
     doc: |2-
        -R --min-supporting-allele-qsum Q Consider any allele in which the sum of qualities of supporting observations is at least Q. default: 0
     type: int
     default: 0
     inputBinding:
       prefix: -R
   - id: minSupMQsum
     label: minSupMQsum
     doc: |2-
        -Y --min-supporting-mapping-qsum Q Consider any allele in which and the sum of mapping qualities of supporting reads is at least Q. default: 0
     type: int
     default: 0
     inputBinding:
       prefix: -Y
   - id: minSupBQthres
     label: minSupBQthres
     doc: |2-
        -Q --mismatch-base-quality-threshold Q Count mismatches toward --read-mismatch-limit if the base quality of the mismatch is >= Q. default: 10
     type: int
     default: 10
     inputBinding:
       prefix: -Q
   - id: readMisMatchLim
     label: readMisMatchLim
     doc: |2-
        -U --read-mismatch-limit N Exclude reads with more than N mismatches where each mismatch has base quality >= mismatch-base-quality-threshold. default: ~unbounded
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -U
   - id: maxMisMatchFrac
     label: maxMisMatchFrac
     doc: |2-
        -z --read-max-mismatch-fraction N Exclude reads with more than N [0,1] fraction of mismatches where each mismatch has base quality >= mismatch-base-quality-threshold default: 1.0
     type: float
     default: 1.0
     inputBinding:
       prefix: -z
   - id: readSNPLim
     label: readSNPLim
     doc: |2-
        -$ --read-snp-limit N Exclude reads with more than N base mismatches, ignoring gaps with quality >= mismatch-base-quality-threshold. default: ~unbounded
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -$
   - id: readINDELLim
     label: readINDELLim
     doc: |2-
        -e --read-indel-limit N Exclude reads with more than N separate gaps. default: ~unbounded
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -e
   - id: standardFilterFlag
     label: standardFilterFlag
     doc: |2-
        -0 --standard-filters Use stringent input base and mapping quality filters Equivalent to -m 30 -q 20 -R 0 -S 0
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-0'
   - id: minAltFrac
     label: minAltFrac
     doc: |2-
        -F --min-alternate-fraction N Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position. default: 0.05
     type: float
     default: 0.05
     inputBinding:
       prefix: -F
   - id: minAltCount
     label: minAltCount
     doc: |2-
        -C --min-alternate-count N Require at least this count of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 2
     type: int
     default: 2
     inputBinding:
       prefix: -C
   - id: minAltQSum
     label: minAltQSum
     doc: |2-
        -3 --min-alternate-qsum N Require at least this sum of quality of observations supporting an alternate allele within a single individual in order to evaluate the position. default: 0
     type: int
     default: 0
     inputBinding:
       prefix: '-3'
   - id: minAltTotal
     label: minAltTotal
     doc: |2-
        -G --min-alternate-total N Require at least this count of observations supporting an alternate allele within the total population in order to use the allele in analysis. default: 1
     type: int
     default: 1
     inputBinding:
       prefix: -G
   - id: minCov
     label: minCov
     doc: ' --min-coverage N Require at least this coverage to process a site. default:
       0'
     type: int
     default: 0
     inputBinding:
       prefix: --min-coverage
   - id: maxCov
     label: maxCov
     doc: |2-
        --max-coverage N Do not process sites with greater than this coverage. default: no limit
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --max-coverage
   - id: noPopPriorsFlag
     label: noPopPriorsFlag
     doc: |2-
        -k --no-population-priors Equivalent to --pooled-discrete --hwe-priors-off and removal of Ewens Sampling Formula component of priors.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -k
   - id: noHWEPriorsFlag
     label: noHWEPriorsFlag
     doc: |2-
        -w --hwe-priors-off Disable estimation of the probability of the combination arising under HWE given the allele frequency as estimated by observation frequency.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -w
   - id: noBinOBSPriorsFlag
     label: noBinOBSPriorsFlag
     doc: |2-
        -V --binomial-obs-priors-off Disable incorporation of prior expectations about observations. Uses read placement probability, strand balance probability, and read position (5'-3') probability.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -V
   - id: noABPriorsFlag
     label: noABPriorsFlag
     doc: |2-
        -a --allele-balance-priors-off Disable use of aggregate probability of observation balance between alleles as a component of the priors.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -a
   - id: obsBiasFile
     label: obsBiasFile
     doc: |2-
        --observation-bias FILE Read length-dependent allele observation biases from FILE. The format is [length] [alignment efficiency relative to reference] where the efficiency is 1 if there is no relative observation bias.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --observation-bias
   - id: baseQualCap
     label: baseQualCap
     doc: |2-
        --base-quality-cap Q Limit estimated observation quality by capping base quality at Q.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --base-quality-cap
   - id: probContamin
     label: probContamin
     doc: |2-
        --prob-contamination F An estimate of contamination to use for all samples. default: 10e-9
     type: float
     default: 1e-09
     inputBinding:
       prefix: --prob-contamination
   - id: legGLScalc
     label: legGLScalc
     doc: ' --legacy-gls Use legacy (polybayes equivalent) genotype likelihood calculations'
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --legacy-gls
   - id: contaminEst
     label: contaminEst
     doc: |2-
        --contamination-estimates FILE A file containing per-sample estimates of contamination, such as those generated by VerifyBamID. The format should be: sample p(read=R|genotype=AR) p(read=A|genotype=AA) Sample '*' can be used to set default contamination estimates.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --contamination-estimates
   - id: repoprtMaxGLFlag
     label: repoprtMaxGLFlag
     doc: |2-
        --report-genotype-likelihood-max Report genotypes using the maximum-likelihood estimate provided from genotype likelihoods.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --report-genotype-likelihood-max
   - id: genotypingMaxIter
     label: genotypingMaxIter
     doc: |2-
        -B --genotyping-max-iterations N Iterate no more than N times during genotyping step. default: 1000.
     type: int
     default: 1000
     inputBinding:
       prefix: -B
   - id: genotypingMaxBDepth
     label: genotypingMaxBDepth
     doc: |2-
        --genotyping-max-banddepth N Integrate no deeper than the Nth best genotype by likelihood when genotyping. default: 6.
     type: int
     default: 6
     inputBinding:
       prefix: --genotyping-max-banddepth
   - id: postIntegrationLim
     label: postIntegrationLim
     doc: |2-
        -W --posterior-integration-limits N,M Integrate all genotype combinations in our posterior space which include no more than N samples with their Mth best data likelihood. default: 1,3.
     type: string
     default: 1,3
     inputBinding:
       prefix: -W
   - id: excludeUnObsGT
     label: excludeUnObsGT
     doc: |2-
        -N --exclude-unobserved-genotypes Skip sample genotypings for which the sample has no supporting reads.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -N
   - id: gtVarThres
     label: gtVarThres
     doc: |2-
        -S --genotype-variant-threshold N Limit posterior integration to samples where the second-best genotype likelihood is no more than log(N) from the highest genotype likelihood for the sample. default: ~unbounded
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -S
   - id: useMQFlag
     label: useMQFlag
     doc: |2-
        -j --use-mapping-quality Use mapping quality of alleles when calculating data likelihoods.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -j
   - id: harmIndelQualFlag
     label: harmIndelQualFlag
     doc: |2-
        -H --harmonic-indel-quality Use a weighted sum of base qualities around an indel, scaled by the distance from the indel. By default use a minimum BQ in flanking sequence.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -H
   - id: readDepFact
     label: readDepFact
     doc: |2-
        -D --read-dependence-factor N Incorporate non-independence of reads by scaling successive observations by this factor during data likelihood calculations. default: 0.9
     type: float
     default: 0.9
     inputBinding:
       prefix: -D
   - id: gtQuals
     label: gtQuals
     doc: |2-
        -= --genotype-qualities Calculate the marginal probability of genotypes and report as GQ in each sample field in the VCF output.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -=

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: freebayes
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: freebayes


