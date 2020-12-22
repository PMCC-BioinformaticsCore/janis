:orphan:

VarDict (Germline)
=====================================

``vardict_germline`` · *1 contributor · 5 versions*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.vardict.vardictgermline import VarDictGermline_1_5_8

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vardict_germline_step",
           VarDictGermline_1_5_8(
               intervals=None,
               bam=None,
               reference=None,
               sampleName=None,
               var2vcfSampleName=None,
               var2vcfAlleleFreqThreshold=None,
           )
       )
       wf.output("out", source=vardict_germline_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vardict_germline:

.. code-block:: bash

   # user inputs
   janis inputs vardict_germline > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       bam: bam.bam
       intervals: intervals.bed
       reference: reference.fasta
       sampleName: <value>
       var2vcfAlleleFreqThreshold: 0.0
       var2vcfSampleName: <value>




5. Run vardict_germline with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vardict_germline





Information
------------

:ID: ``vardict_germline``
:URL: *No URL to the documentation was provided*
:Versions: 1.7.0, 1.6.0, 1.5.8, 1.5.7, 1.5.6
:Container: michaelfranklin/vardict:1.5.8
:Authors: Michael Franklin
:Citations: None
:Created: 2019-01-21
:Updated: 2020-06-04


Outputs
-----------

======  ======  ===============
name    type    documentation
======  ======  ===============
out     VCF
======  ======  ===============


Additional configuration (inputs)
---------------------------------

==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================
name                        type                prefix      position  documentation
==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================
intervals                   bed                                    2
bam                         IndexedBam          -b                 1  The indexed BAM file
reference                   FastaFai            -G                 1  The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
sampleName                  String              -N                 1  The sample name to be used directly.  Will overwrite -n option
var2vcfSampleName           String              -N                 5
var2vcfAlleleFreqThreshold  Float               -f                 5
outputFilename              Optional<Filename>  >                  6
indels3prime                Optional<Boolean>   -3                 1  Indicate to move indels to 3-prime if alternative alignment can be achieved.
amplicon                    Optional<Float>     -a                 1  Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
minReads                    Optional<Integer>   -B                 1  The minimum # of reads to determine strand bias, default 2
chromNamesAreNumbers        Optional<Boolean>   -C                 1  Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
chromColumn                 Optional<Integer>   -c                 1  The column for chromosome
debug                       Optional<Boolean>   -D                 1  Debug mode.  Will print some error messages and append full genotype at the end.
splitDelimeter              Optional<String>    -d                 1  The delimiter for split region_info, default to tab "	"
geneEndCol                  Optional<Integer>   -E                 1  The column for region end, e.g. gene end
segEndCol                   Optional<Integer>   -e                 1  The column for segment ends in the region, e.g. exon ends
filter                      Optional<String>    -F                 1  The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
alleleFreqThreshold         Optional<Float>     -f                 1  The threshold for allele frequency, default: 0.05 or 5%
geneNameCol                 Optional<Integer>   -g                 1  The column for gene name, or segment annotation
printHeaderRow              Optional<Boolean>   -h                 1  Print a header row describing columns
indelSize                   Optional<Integer>   -I                 1  The indel size.  Default: 120bp
outputSplice                Optional<Boolean>   -i                 1  Output splicing read counts
performLocalRealignment     Optional<Integer>   -k                 1  Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
minMatches                  Optional<Integer>   -M                 1  The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
maxMismatches               Optional<Integer>   -m                 1  If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
regexSampleName             Optional<String>    -n                 1  The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
mapq                        Optional<String>    -O                 1  The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
qratio                      Optional<Float>     -o                 1  The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
readPosition                Optional<Float>     -P                 1  The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
pileup                      Optional<Boolean>   -p                 1  Do pileup regardless of the frequency
minMappingQual              Optional<Integer>   -Q                 1  If set, reads with mapping quality less than INT will be filtered and ignored
phredScore                  Optional<Integer>   -q                 1  The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
region                      Optional<String>    -R                 1  The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
minVariantReads             Optional<Integer>   -r                 1  The minimum # of variant reads, default 2
regStartCol                 Optional<Integer>   -S                 1  The column for region start, e.g. gene start
segStartCol                 Optional<Integer>   -s                 1  The column for segment starts in the region, e.g. exon starts
minReadsBeforeTrim          Optional<Integer>   -T                 1  Trim bases after [INT] bases in the reads
removeDuplicateReads        Optional<Boolean>   -t                 1  Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
threads                     Optional<Integer>   -th                1  Threads count.
freq                        Optional<Integer>   -V                 1  The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
vcfFormat                   Optional<Boolean>   -v                 1  VCF format output
vs                          Optional<String>    -VS                1  [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
bp                          Optional<Integer>   -X                 1  Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
extensionNucleotide         Optional<Integer>   -x                 1  The number of nucleotide to extend for each segment, default: 0
yy                          Optional<Boolean>   -y                 1  <No content>
downsamplingFraction        Optional<Integer>   -Z                 1  For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
zeroBasedCoords             Optional<Integer>   -z                 1  0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
==========================  ==================  ========  ==========  ==================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vardict_germline {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File intervals
       String? outputFilename
       File bam
       File bam_bai
       File reference
       File reference_fai
       Boolean? indels3prime
       Float? amplicon
       Int? minReads
       Boolean? chromNamesAreNumbers
       Int? chromColumn
       Boolean? debug
       String? splitDelimeter
       Int? geneEndCol
       Int? segEndCol
       String? filter
       Float? alleleFreqThreshold
       Int? geneNameCol
       Boolean? printHeaderRow
       Int? indelSize
       Boolean? outputSplice
       Int? performLocalRealignment
       Int? minMatches
       Int? maxMismatches
       String sampleName
       String? regexSampleName
       String? mapq
       Float? qratio
       Float? readPosition
       Boolean? pileup
       Int? minMappingQual
       Int? phredScore
       String? region
       Int? minVariantReads
       Int? regStartCol
       Int? segStartCol
       Int? minReadsBeforeTrim
       Boolean? removeDuplicateReads
       Int? threads
       Int? freq
       Boolean? vcfFormat
       String? vs
       Int? bp
       Int? extensionNucleotide
       Boolean? yy
       Int? downsamplingFraction
       Int? zeroBasedCoords
       String var2vcfSampleName
       Float var2vcfAlleleFreqThreshold
     }
     command <<<
       set -e
       VarDict \
         -b ~{bam} \
         -G ~{reference} \
         ~{if (defined(indels3prime) && select_first([indels3prime])) then "-3" else ""} \
         ~{if defined(amplicon) then ("-a " + amplicon) else ''} \
         ~{if defined(minReads) then ("-B " + minReads) else ''} \
         ~{if (defined(chromNamesAreNumbers) && select_first([chromNamesAreNumbers])) then "-C" else ""} \
         ~{if defined(chromColumn) then ("-c " + chromColumn) else ''} \
         ~{if (defined(debug) && select_first([debug])) then "-D" else ""} \
         ~{if defined(splitDelimeter) then ("-d " + splitDelimeter) else ''} \
         ~{if defined(geneEndCol) then ("-E " + geneEndCol) else ''} \
         ~{if defined(segEndCol) then ("-e " + segEndCol) else ''} \
         ~{if defined(filter) then ("-F " + filter) else ''} \
         ~{if defined(alleleFreqThreshold) then ("-f " + alleleFreqThreshold) else ''} \
         ~{if defined(geneNameCol) then ("-g " + geneNameCol) else ''} \
         ~{if (defined(printHeaderRow) && select_first([printHeaderRow])) then "-h" else ""} \
         ~{if defined(indelSize) then ("-I " + indelSize) else ''} \
         ~{if (defined(outputSplice) && select_first([outputSplice])) then "-i" else ""} \
         ~{if defined(performLocalRealignment) then ("-k " + performLocalRealignment) else ''} \
         ~{if defined(minMatches) then ("-M " + minMatches) else ''} \
         ~{if defined(maxMismatches) then ("-m " + maxMismatches) else ''} \
         -N ~{sampleName} \
         ~{if defined(regexSampleName) then ("-n " + regexSampleName) else ''} \
         ~{if defined(mapq) then ("-O " + mapq) else ''} \
         ~{if defined(qratio) then ("-o " + qratio) else ''} \
         ~{if defined(readPosition) then ("-P " + readPosition) else ''} \
         ~{if (defined(pileup) && select_first([pileup])) then "-p" else ""} \
         ~{if defined(minMappingQual) then ("-Q " + minMappingQual) else ''} \
         ~{if defined(phredScore) then ("-q " + phredScore) else ''} \
         ~{if defined(region) then ("-R " + region) else ''} \
         ~{if defined(minVariantReads) then ("-r " + minVariantReads) else ''} \
         ~{if defined(regStartCol) then ("-S " + regStartCol) else ''} \
         ~{if defined(segStartCol) then ("-s " + segStartCol) else ''} \
         ~{if defined(minReadsBeforeTrim) then ("-T " + minReadsBeforeTrim) else ''} \
         ~{if (defined(removeDuplicateReads) && select_first([removeDuplicateReads])) then "-t" else ""} \
         ~{if defined(select_first([threads, select_first([runtime_cpu, 1])])) then ("-th " + select_first([threads, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(freq) then ("-V " + freq) else ''} \
         ~{if (defined(vcfFormat) && select_first([vcfFormat])) then "-v" else ""} \
         ~{if defined(vs) then ("-VS " + vs) else ''} \
         ~{if defined(bp) then ("-X " + bp) else ''} \
         ~{if defined(extensionNucleotide) then ("-x " + extensionNucleotide) else ''} \
         ~{if (defined(yy) && select_first([yy])) then "-y" else ""} \
         ~{if defined(downsamplingFraction) then ("-Z " + downsamplingFraction) else ''} \
         ~{if defined(zeroBasedCoords) then ("-z " + zeroBasedCoords) else ''} \
         ~{intervals} \
         | teststrandbias.R | \
         var2vcf_valid.pl \
         -N ~{var2vcfSampleName} \
         -f ~{var2vcfAlleleFreqThreshold} \
         > ~{select_first([outputFilename, "generated.vardict.vcf"])}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "michaelfranklin/vardict:1.5.8"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 8, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.vardict.vcf"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: VarDict (Germline)
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: michaelfranklin/vardict:1.5.8

   inputs:
   - id: intervals
     label: intervals
     type: File
     inputBinding:
       position: 2
       shellQuote: false
   - id: outputFilename
     label: outputFilename
     type:
     - string
     - 'null'
     default: generated.vardict.vcf
     inputBinding:
       prefix: '>'
       position: 6
       shellQuote: false
   - id: bam
     label: bam
     doc: The indexed BAM file
     type: File
     secondaryFiles:
     - pattern: .bai
     inputBinding:
       prefix: -b
       position: 1
       shellQuote: false
   - id: reference
     label: reference
     doc: |-
       The reference fasta. Should be indexed (.fai). Defaults to: /ngs/reference_data/genomes/Hsapiens/hg19/seq/hg19.fa
     type: File
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: -G
       position: 1
       shellQuote: false
   - id: indels3prime
     label: indels3prime
     doc: Indicate to move indels to 3-prime if alternative alignment can be achieved.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: '-3'
       position: 1
       shellQuote: false
   - id: amplicon
     label: amplicon
     doc: |-
       Indicate it's amplicon based calling.  Reads that don't map to the amplicon will be skipped.  A read pair is considered belonging  to the amplicon if the edges are less than int bp to the amplicon, and overlap fraction is at least float.  Default: 10:0.95
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -a
       position: 1
       shellQuote: false
   - id: minReads
     label: minReads
     doc: 'The minimum # of reads to determine strand bias, default 2'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -B
       position: 1
       shellQuote: false
   - id: chromNamesAreNumbers
     label: chromNamesAreNumbers
     doc: Indicate the chromosome names are just numbers, such as 1, 2, not chr1, chr2
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -C
       position: 1
       shellQuote: false
   - id: chromColumn
     label: chromColumn
     doc: The column for chromosome
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -c
       position: 1
       shellQuote: false
   - id: debug
     label: debug
     doc: Debug mode.  Will print some error messages and append full genotype at the
       end.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -D
       position: 1
       shellQuote: false
   - id: splitDelimeter
     label: splitDelimeter
     doc: "The delimiter for split region_info, default to tab \"\t\""
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -d
       position: 1
       shellQuote: false
   - id: geneEndCol
     label: geneEndCol
     doc: The column for region end, e.g. gene end
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -E
       position: 1
       shellQuote: false
   - id: segEndCol
     label: segEndCol
     doc: The column for segment ends in the region, e.g. exon ends
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -e
       position: 1
       shellQuote: false
   - id: filter
     label: filter
     doc: |-
       The hexical to filter reads using samtools. Default: 0x500 (filter 2nd alignments and duplicates). Use -F 0 to turn it off.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -F
       position: 1
       shellQuote: false
   - id: alleleFreqThreshold
     label: alleleFreqThreshold
     doc: 'The threshold for allele frequency, default: 0.05 or 5%'
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -f
       position: 1
       shellQuote: false
   - id: geneNameCol
     label: geneNameCol
     doc: The column for gene name, or segment annotation
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -g
       position: 1
       shellQuote: false
   - id: printHeaderRow
     label: printHeaderRow
     doc: Print a header row describing columns
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -h
       position: 1
       shellQuote: false
   - id: indelSize
     label: indelSize
     doc: 'The indel size.  Default: 120bp'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -I
       position: 1
       shellQuote: false
   - id: outputSplice
     label: outputSplice
     doc: Output splicing read counts
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -i
       position: 1
       shellQuote: false
   - id: performLocalRealignment
     label: performLocalRealignment
     doc: |-
       Indicate whether to perform local realignment.  Default: 1.  Set to 0 to disable it. For Ion or PacBio, 0 is recommended.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -k
       position: 1
       shellQuote: false
   - id: minMatches
     label: minMatches
     doc: |-
       The minimum matches for a read to be considered. If, after soft-clipping, the matched bp is less than INT, then the read is discarded. It's meant for PCR based targeted sequencing where there's no insert and the matching is only the primers. Default: 0, or no filtering
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -M
       position: 1
       shellQuote: false
   - id: maxMismatches
     label: maxMismatches
     doc: |-
       If set, reads with mismatches more than INT will be filtered and ignored. Gaps are not counted as mismatches. Valid only for bowtie2/TopHat or BWA aln followed by sampe. BWA mem is calculated as NM - Indels. Default: 8, or reads with more than 8 mismatches will not be used.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -m
       position: 1
       shellQuote: false
   - id: sampleName
     label: sampleName
     doc: The sample name to be used directly.  Will overwrite -n option
     type: string
     inputBinding:
       prefix: -N
       position: 1
       shellQuote: false
   - id: regexSampleName
     label: regexSampleName
     doc: |-
       The regular expression to extract sample name from BAM filenames. Default to: /([^\/\._]+?)_[^\/]*.bam/
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -n
       position: 1
       shellQuote: false
   - id: mapq
     label: mapq
     doc: |-
       The reads should have at least mean MapQ to be considered a valid variant. Default: no filtering
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -O
       position: 1
       shellQuote: false
   - id: qratio
     label: qratio
     doc: |-
       The Qratio of (good_quality_reads)/(bad_quality_reads+0.5). The quality is defined by -q option.  Default: 1.5
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -o
       position: 1
       shellQuote: false
   - id: readPosition
     label: readPosition
     doc: |-
       The read position filter. If the mean variants position is less that specified, it's considered false positive.  Default: 5
     type:
     - float
     - 'null'
     inputBinding:
       prefix: -P
       position: 1
       shellQuote: false
   - id: pileup
     label: pileup
     doc: Do pileup regardless of the frequency
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -p
       position: 1
       shellQuote: false
   - id: minMappingQual
     label: minMappingQual
     doc: If set, reads with mapping quality less than INT will be filtered and ignored
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -Q
       position: 1
       shellQuote: false
   - id: phredScore
     label: phredScore
     doc: |-
       The phred score for a base to be considered a good call.  Default: 25 (for Illumina) For PGM, set it to ~15, as PGM tends to under estimate base quality.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -q
       position: 1
       shellQuote: false
   - id: region
     label: region
     doc: |-
       The region of interest.  In the format of chr:start-end.  If end is omitted, then a single position.  No BED is needed.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -R
       position: 1
       shellQuote: false
   - id: minVariantReads
     label: minVariantReads
     doc: 'The minimum # of variant reads, default 2'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -r
       position: 1
       shellQuote: false
   - id: regStartCol
     label: regStartCol
     doc: The column for region start, e.g. gene start
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -S
       position: 1
       shellQuote: false
   - id: segStartCol
     label: segStartCol
     doc: The column for segment starts in the region, e.g. exon starts
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -s
       position: 1
       shellQuote: false
   - id: minReadsBeforeTrim
     label: minReadsBeforeTrim
     doc: Trim bases after [INT] bases in the reads
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -T
       position: 1
       shellQuote: false
   - id: removeDuplicateReads
     label: removeDuplicateReads
     doc: |-
       Indicate to remove duplicated reads.  Only one pair with same start positions will be kept
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -t
       position: 1
       shellQuote: false
   - id: threads
     label: threads
     doc: Threads count.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -th
       position: 1
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
       shellQuote: false
   - id: freq
     label: freq
     doc: |-
       The lowest frequency in the normal sample allowed for a putative somatic mutation. Defaults to 0.05
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -V
       position: 1
       shellQuote: false
   - id: vcfFormat
     label: vcfFormat
     doc: VCF format output
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -v
       position: 1
       shellQuote: false
   - id: vs
     label: vs
     doc: |-
       [STRICT | LENIENT | SILENT] How strict to be when reading a SAM or BAM: STRICT   - throw an exception if something looks wrong. LENIENT	- Emit warnings but keep going if possible. SILENT	- Like LENIENT, only don't emit warning messages. Default: LENIENT
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -VS
       position: 1
       shellQuote: false
   - id: bp
     label: bp
     doc: |-
       Extension of bp to look for mismatches after insersion or deletion.  Default to 3 bp, or only calls when they're within 3 bp.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -X
       position: 1
       shellQuote: false
   - id: extensionNucleotide
     label: extensionNucleotide
     doc: 'The number of nucleotide to extend for each segment, default: 0'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -x
       position: 1
       shellQuote: false
   - id: yy
     label: yy
     doc: <No content>
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -y
       position: 1
       shellQuote: false
   - id: downsamplingFraction
     label: downsamplingFraction
     doc: |-
       For downsampling fraction.  e.g. 0.7 means roughly 70% downsampling.  Default: No downsampling.  Use with caution.  The downsampling will be random and non-reproducible.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -Z
       position: 1
       shellQuote: false
   - id: zeroBasedCoords
     label: zeroBasedCoords
     doc: |-
       0/1  Indicate whether coordinates are zero-based, as IGV uses.  Default: 1 for BED file or amplicon BED file. Use 0 to turn it off. When using the -R option, it's set to 0
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -z
       position: 1
       shellQuote: false
   - id: var2vcfSampleName
     label: var2vcfSampleName
     type: string
     inputBinding:
       prefix: -N
       position: 5
       shellQuote: false
   - id: var2vcfAlleleFreqThreshold
     label: var2vcfAlleleFreqThreshold
     type: float
     inputBinding:
       prefix: -f
       position: 5
       shellQuote: false

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.vardict.vcf
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: VarDict
   arguments:
   - position: 3
     valueFrom: '| teststrandbias.R |'
     shellQuote: false
   - position: 4
     valueFrom: var2vcf_valid.pl
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vardict_germline


