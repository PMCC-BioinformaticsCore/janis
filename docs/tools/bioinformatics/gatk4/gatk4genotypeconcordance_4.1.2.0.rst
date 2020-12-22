:orphan:

GATK4: Genotype Concordance
======================================================

``Gatk4GenotypeConcordance`` · *1 contributor · 4 versions*

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


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.gatk4.genotypeconcordance.versions import Gatk4GenotypeConcordance_4_1_2

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "gatk4genotypeconcordance_step",
           Gatk4GenotypeConcordance_4_1_2(
               callVCF=None,
               truthVCF=None,
           )
       )
       wf.output("summaryMetrics", source=gatk4genotypeconcordance_step.summaryMetrics)
       wf.output("detailMetrics", source=gatk4genotypeconcordance_step.detailMetrics)
       wf.output("contingencyMetrics", source=gatk4genotypeconcordance_step.contingencyMetrics)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for Gatk4GenotypeConcordance:

.. code-block:: bash

   # user inputs
   janis inputs Gatk4GenotypeConcordance > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       callVCF: callVCF.vcf.gz
       truthVCF: truthVCF.vcf




5. Run Gatk4GenotypeConcordance with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       Gatk4GenotypeConcordance





Information
------------

:ID: ``Gatk4GenotypeConcordance``
:URL: `https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php <https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.5.0/picard_vcf_GenotypeConcordance.php>`_
:Versions: 4.1.4.0, 4.1.3.0, 4.1.2.0, 4.0.12.0
:Container: broadinstitute/gatk:4.1.2.0
:Authors: Michael Franklin
:Citations: See https://software.broadinstitute.org/gatk/documentation/article?id=11027 for more information
:Created: 2018-12-24
:Updated: 2019-01-24


Outputs
-----------

==================  ======  ===============
name                type    documentation
==================  ======  ===============
summaryMetrics      File
detailMetrics       File
contingencyMetrics  File
==================  ======  ===============


Additional configuration (inputs)
---------------------------------

==========================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
name                        type                     prefix                     position  documentation
==========================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================
callVCF                     Gzipped<VCF>             --CALL_VCF                           The VCF containing the call sample
truthVCF                    IndexedVCF               --TRUTH_VCF                          The VCF containing the truth sample
javaOptions                 Optional<Array<String>>
compression_level           Optional<Integer>                                             Compression level for all compressed files created (e.g. BAM and VCF). Default value: 2.
outputBasename              Optional<Filename>       --OUTPUT                             Basename for the three metrics files that are to be written. Resulting files will be:(1) .genotype_concordance_summary_metrics, (2) .genotype_concordance_detail_metrics, (3) .genotype_concordance_contingency_metrics.
argumentsFile               Optional<Array<File>>    --arguments_file                 10  read one or more arguments files and add them to the command line
callSample                  Optional<String>         --CALL_SAMPLE                    10  The name of the call sample within the call VCF. Not required if only one sample exists.
ignoreFilterStatus          Optional<Boolean>        --IGNORE_FILTER_STATUS               Default is false. If true, filter status of sites will be ignored so that we include filtered sites when calculating genotype concordance.
intersectIntervals          Optional<Boolean>        --INTERSECT_INTERVALS                If true, multiple interval lists will be intersected. If false multiple lists will be unioned.
intervals                   Optional<Array<VCF>>     --INTERVALS                          One or more interval list files that will be used to limit the genotype concordance. Note - if intervals are specified, the VCF files must be indexed.
minDP                       Optional<Float>          --MIN_DP                             Genotypes below this depth will have genotypes classified as LowDp.
minGQ                       Optional<Float>          --MIN_GQ                             Genotypes below this genotype quality will have genotypes classified as LowGq.
treatMissingSitesAsHomeRef  Optional<Boolean>        --MISSING_SITES_HOM_REF              Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth
                                                                                          set will be treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. This flag can only be used with a high confidence interval list.
outputAllRows               Optional<Boolean>        --OUTPUT_ALL_ROWS                    If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.
outputVcf                   Optional<Boolean>        --OUTPUT_VCF                         Output a VCF annotated with concordance information.
truthSample                 Optional<String>         --TRUTH_SAMPLE                       The name of the truth sample within the truth VCF. Not required if only one sample exists.
useVcfIndex                 Optional<Boolean>        --USE_VCF_INDEX                      If true, use the VCF index, else iterate over the entire VCF
compressionLevel            Optional<Integer>        --COMPRESSION_LEVEL              11  Compression level for all compressed files created (e.g. BAM and GELI).
createIndex                 Optional<Boolean>        --CREATE_INDEX                   11  Whether to create a BAM index when writing a coordinate-sorted BAM file.
createMd5File               Optional<Boolean>        --CREATE_MD5_FILE                11  Whether to create an MD5 digest for any BAM or FASTQ files created.
maxRecordsInRam             Optional<Integer>        --MAX_RECORDS_IN_RAM             11  When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
quiet                       Optional<Boolean>        --QUIET                          11  Whether to suppress job-summary info on System.err.
reference                   Optional<File>           --REFERENCE_SEQUENCE             11  Reference sequence file.
tmpDir                      Optional<String>         --TMP_DIR                        11  Undocumented option
useJdkDeflater              Optional<Boolean>        --use_jdk_deflater               11  Whether to use the JdkDeflater (as opposed to IntelDeflater)
useJdkInflater              Optional<Boolean>        --use_jdk_inflater               11  Whether to use the JdkInflater (as opposed to IntelInflater)
validationStringency        Optional<String>         --VALIDATION_STRINGENCY          11  Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
verbosity                   Optional<String>         --verbosity                      11  The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
==========================  =======================  =======================  ==========  ================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task Gatk4GenotypeConcordance {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       Array[String]? javaOptions
       Int? compression_level
       File callVCF
       File callVCF_tbi
       File truthVCF
       File truthVCF_idx
       String? outputBasename
       Array[File]? argumentsFile
       String? callSample
       Boolean? ignoreFilterStatus
       Boolean? intersectIntervals
       Array[File]? intervals
       Float? minDP
       Float? minGQ
       Boolean? treatMissingSitesAsHomeRef
       Boolean? outputAllRows
       Boolean? outputVcf
       String? truthSample
       Boolean? useVcfIndex
       Int? compressionLevel
       Boolean? createIndex
       Boolean? createMd5File
       Int? maxRecordsInRam
       Boolean? quiet
       File? reference
       String? tmpDir
       Boolean? useJdkDeflater
       Boolean? useJdkInflater
       String? validationStringency
       String? verbosity
     }
     command <<<
       set -e
       gatk GenotypeConcordance \
         --java-options '-Xmx~{((select_first([runtime_memory, 4]) * 3) / 4)}G ~{if (defined(compression_level)) then ("-Dsamjdk.compress_level=" + compression_level) else ""} ~{sep(" ", select_first([javaOptions, []]))}' \
         --CALL_VCF '~{callVCF}' \
         --TRUTH_VCF '~{truthVCF}' \
         --OUTPUT '~{select_first([outputBasename, "generated"])}' \
         ~{if (defined(ignoreFilterStatus) && select_first([ignoreFilterStatus])) then "--IGNORE_FILTER_STATUS" else ""} \
         ~{if (defined(intersectIntervals) && select_first([intersectIntervals])) then "--INTERSECT_INTERVALS" else ""} \
         ~{if (defined(intervals) && length(select_first([intervals])) > 0) then "--INTERVALS '" + sep("' '", select_first([intervals])) + "'" else ""} \
         ~{if defined(minDP) then ("--MIN_DP " + minDP) else ''} \
         ~{if defined(minGQ) then ("--MIN_GQ " + minGQ) else ''} \
         ~{if (defined(treatMissingSitesAsHomeRef) && select_first([treatMissingSitesAsHomeRef])) then "--MISSING_SITES_HOM_REF" else ""} \
         ~{if (defined(outputAllRows) && select_first([outputAllRows])) then "--OUTPUT_ALL_ROWS" else ""} \
         ~{if (defined(outputVcf) && select_first([outputVcf])) then "--OUTPUT_VCF" else ""} \
         ~{if defined(truthSample) then ("--TRUTH_SAMPLE '" + truthSample + "'") else ""} \
         ~{if (defined(useVcfIndex) && select_first([useVcfIndex])) then "--USE_VCF_INDEX" else ""} \
         ~{if (defined(argumentsFile) && length(select_first([argumentsFile])) > 0) then "--arguments_file '" + sep("' '", select_first([argumentsFile])) + "'" else ""} \
         ~{if defined(callSample) then ("--CALL_SAMPLE '" + callSample + "'") else ""} \
         ~{if defined(compressionLevel) then ("--COMPRESSION_LEVEL " + compressionLevel) else ''} \
         ~{if (defined(createIndex) && select_first([createIndex])) then "--CREATE_INDEX" else ""} \
         ~{if (defined(createMd5File) && select_first([createMd5File])) then "--CREATE_MD5_FILE" else ""} \
         ~{if defined(maxRecordsInRam) then ("--MAX_RECORDS_IN_RAM " + maxRecordsInRam) else ''} \
         ~{if (defined(quiet) && select_first([quiet])) then "--QUIET" else ""} \
         ~{if defined(reference) then ("--REFERENCE_SEQUENCE '" + reference + "'") else ""} \
         ~{if defined(select_first([tmpDir, "/tmp/"])) then ("--TMP_DIR '" + select_first([tmpDir, "/tmp/"]) + "'") else ""} \
         ~{if (defined(useJdkDeflater) && select_first([useJdkDeflater])) then "--use_jdk_deflater" else ""} \
         ~{if (defined(useJdkInflater) && select_first([useJdkInflater])) then "--use_jdk_inflater" else ""} \
         ~{if defined(validationStringency) then ("--VALIDATION_STRINGENCY '" + validationStringency + "'") else ""} \
         ~{if defined(verbosity) then ("--verbosity '" + verbosity + "'") else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "broadinstitute/gatk:4.1.2.0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File summaryMetrics = glob("*.genotype_concordance_summary_metrics")[0]
       File detailMetrics = glob("*.genotype_concordance_detail_metrics")[0]
       File contingencyMetrics = glob("*.genotype_concordance_contingency_metrics")[0]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'GATK4: Genotype Concordance'
   doc: |-
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

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: broadinstitute/gatk:4.1.2.0

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
   - id: callVCF
     label: callVCF
     doc: The VCF containing the call sample
     type: File
     secondaryFiles:
     - pattern: .tbi
     inputBinding:
       prefix: --CALL_VCF
   - id: truthVCF
     label: truthVCF
     doc: The VCF containing the truth sample
     type: File
     secondaryFiles:
     - pattern: .idx
     inputBinding:
       prefix: --TRUTH_VCF
   - id: outputBasename
     label: outputBasename
     doc: |-
       Basename for the three metrics files that are to be written. Resulting files will be:(1) .genotype_concordance_summary_metrics, (2) .genotype_concordance_detail_metrics, (3) .genotype_concordance_contingency_metrics.
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --OUTPUT
   - id: argumentsFile
     label: argumentsFile
     doc: read one or more arguments files and add them to the command line
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --arguments_file
       position: 10
   - id: callSample
     label: callSample
     doc: |-
       The name of the call sample within the call VCF. Not required if only one sample exists.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --CALL_SAMPLE
       position: 10
   - id: ignoreFilterStatus
     label: ignoreFilterStatus
     doc: |-
       Default is false. If true, filter status of sites will be ignored so that we include filtered sites when calculating genotype concordance.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --IGNORE_FILTER_STATUS
   - id: intersectIntervals
     label: intersectIntervals
     doc: |-
       If true, multiple interval lists will be intersected. If false multiple lists will be unioned.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --INTERSECT_INTERVALS
   - id: intervals
     label: intervals
     doc: |-
       One or more interval list files that will be used to limit the genotype concordance. Note - if intervals are specified, the VCF files must be indexed.
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --INTERVALS
   - id: minDP
     label: minDP
     doc: Genotypes below this depth will have genotypes classified as LowDp.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --MIN_DP
   - id: minGQ
     label: minGQ
     doc: Genotypes below this genotype quality will have genotypes classified as LowGq.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --MIN_GQ
   - id: treatMissingSitesAsHomeRef
     label: treatMissingSitesAsHomeRef
     doc: |-
       Default is false, which follows the GA4GH Scheme. If true, missing sites in the truth 
       set will be treated as HOM_REF sites and sites missing in both the truth and call sets will be true negatives. Useful when hom ref sites are left out of the truth set. This flag can only be used with a high confidence interval list.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --MISSING_SITES_HOM_REF
   - id: outputAllRows
     label: outputAllRows
     doc: |-
       If true, output all rows in detailed statistics even when count == 0. When false only output rows with non-zero counts.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --OUTPUT_ALL_ROWS
   - id: outputVcf
     label: outputVcf
     doc: Output a VCF annotated with concordance information.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --OUTPUT_VCF
   - id: truthSample
     label: truthSample
     doc: |-
       The name of the truth sample within the truth VCF. Not required if only one sample exists.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --TRUTH_SAMPLE
   - id: useVcfIndex
     label: useVcfIndex
     doc: If true, use the VCF index, else iterate over the entire VCF
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --USE_VCF_INDEX
   - id: compressionLevel
     label: compressionLevel
     doc: Compression level for all compressed files created (e.g. BAM and GELI).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --COMPRESSION_LEVEL
       position: 11
   - id: createIndex
     label: createIndex
     doc: Whether to create a BAM index when writing a coordinate-sorted BAM file.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_INDEX
       position: 11
   - id: createMd5File
     label: createMd5File
     doc: Whether to create an MD5 digest for any BAM or FASTQ files created.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --CREATE_MD5_FILE
       position: 11
   - id: maxRecordsInRam
     label: maxRecordsInRam
     doc: |-
       When writing SAM files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort a SAM file, and increases the amount of RAM needed.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --MAX_RECORDS_IN_RAM
       position: 11
   - id: quiet
     label: quiet
     doc: Whether to suppress job-summary info on System.err.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --QUIET
       position: 11
   - id: reference
     label: reference
     doc: Reference sequence file.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --REFERENCE_SEQUENCE
       position: 11
   - id: tmpDir
     label: tmpDir
     doc: Undocumented option
     type: string
     default: /tmp/
     inputBinding:
       prefix: --TMP_DIR
       position: 11
   - id: useJdkDeflater
     label: useJdkDeflater
     doc: Whether to use the JdkDeflater (as opposed to IntelDeflater)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_jdk_deflater
       position: 11
   - id: useJdkInflater
     label: useJdkInflater
     doc: Whether to use the JdkInflater (as opposed to IntelInflater)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_jdk_inflater
       position: 11
   - id: validationStringency
     label: validationStringency
     doc: |-
       Validation stringency for all SAM files read by this program. Setting stringency to SILENT can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) do not otherwise need to be decoded.The --VALIDATION_STRINGENCY argument is an enumerated type (ValidationStringency), which can have one of the following values: [STRICT, LENIENT, SILENT]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --VALIDATION_STRINGENCY
       position: 11
   - id: verbosity
     label: verbosity
     doc: |-
       The --verbosity argument is an enumerated type (LogLevel), which can have one of the following values: [ERROR, WARNING, INFO, DEBUG]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --verbosity
       position: 11

   outputs:
   - id: summaryMetrics
     label: summaryMetrics
     type: File
     outputBinding:
       glob: '*.genotype_concordance_summary_metrics'
       loadContents: false
   - id: detailMetrics
     label: detailMetrics
     type: File
     outputBinding:
       glob: '*.genotype_concordance_detail_metrics'
       loadContents: false
   - id: contingencyMetrics
     label: contingencyMetrics
     type: File
     outputBinding:
       glob: '*.genotype_concordance_contingency_metrics'
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - gatk
   - GenotypeConcordance
   arguments:
   - prefix: --java-options
     position: -1
     valueFrom: |-
       $("-Xmx{memory}G {compression} {otherargs}".replace(/\{memory\}/g, (([inputs.runtime_memory, 4].filter(function (inner) { return inner != null })[0] * 3) / 4)).replace(/\{compression\}/g, (inputs.compression_level != null) ? ("-Dsamjdk.compress_level=" + inputs.compression_level) : "").replace(/\{otherargs\}/g, [inputs.javaOptions, []].filter(function (inner) { return inner != null })[0].join(" ")))

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: Gatk4GenotypeConcordance


