:orphan:

BCFTools: Normalize
==================================

``bcftoolsNorm`` · *1 contributor · 2 versions*

Left-align and normalize indels, check if REF alleles match the reference,
split multiallelic sites into multiple rows; recover multiallelics from
multiple rows. Left-alignment and normalization will only be applied if
the --fasta-ref option is supplied.



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.bcftools.norm.versions import BcfToolsNorm_1_5

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "bcftoolsnorm_step",
           BcfToolsNorm_1_5(
               vcf=None,
           )
       )
       wf.output("out", source=bcftoolsnorm_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for bcftoolsNorm:

.. code-block:: bash

   # user inputs
   janis inputs bcftoolsNorm > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       vcf: vcf.vcf.gz




5. Run bcftoolsNorm with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       bcftoolsNorm





Information
------------

:ID: ``bcftoolsNorm``
:URL: `https://samtools.github.io/bcftools/bcftools.html#norm <https://samtools.github.io/bcftools/bcftools.html#norm>`_
:Versions: v1.9, v1.5
:Container: biocontainers/bcftools:v1.5_cv2
:Authors: Michael Franklin
:Citations: Li H, Handsaker B, Wysoker A, Fennell T, Ruan J, Homer N, Marth G, Abecasis G, Durbin R, and 1000 Genome Project Data Processing Subgroup, The Sequence alignment/map (SAM) format and SAMtools, Bioinformatics (2009) 25(16) 2078-9
:DOI: http://www.ncbi.nlm.nih.gov/pubmed/19505943
:Created: 2019-01-24
:Updated: 2019-01-24


Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
out     Gzipped<VCF>
======  ============  ===============


Additional configuration (inputs)
---------------------------------

=====================  =====================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                   type                   prefix          position  documentation
=====================  =====================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
vcf                    Gzipped<VCF>                                 10
outputFilename         Optional<Filename>     -o                        --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
checkRef               Optional<String>       -c                        --check-ref e|w|x|s: what to do when incorrect or missing REF allele is encountered: exit (e), warn (w), exclude (x), or set/fix (s) bad sites. The w option can be combined with x and s. Note that s can swap alleles and will update genotypes (GT) and AC counts, but will not attempt to fix PL or other fields. Also note, and this cannot be stressed enough, that s will NOT fix strand issues in your VCF, do NOT use it for that purpose!!! (Instead see http://samtools.github.io/bcftools/howtos/plugin.af-dist.html and http://samtools.github.io/bcftools/howtos/plugin.fixref.html.)
removeDups             Optional<String>       -d                        --rm-dup: snps|indels|both|all|none. If a record is present multiple times, output only the first instance, see --collapse in Common Options.
removeDupsAcrossFiles  Optional<Boolean>      -D                        --remove-duplicates: If a record is present in multiple files, output only the first instance. Alias for -d none, deprecated.
reference              Optional<FastaFai>     -f                        --fasta-ref: reference sequence. Supplying this option will turn on left-alignment and normalization, however, see also the --do-not-normalize option below.
multiallelics          Optional<String>       -m                        --multiallelics -|+[snps|indels|both|any]: split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
noVersion              Optional<Boolean>      --no-version              Do not append version and command line information to the output VCF header.
noNormalize            Optional<Boolean>      -N                        --do-not-normalize: the -c s option can be used to fix or set the REF allele from the reference -f. The -N option will not turn on indel normalisation as the -f option normally implies
outputType             Optional<String>       -O                        --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
regions                Optional<String>       -r                        --regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.
regionsFile            Optional<File>         -R                        --regions-file: Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20' is not the same as '20'. Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.
strictFilter           Optional<Boolean>      -s                        --strict-filter: when merging (-m+), merged site is PASS only if all sites being merged PASS
targets                Optional<Array<File>>  -t                        --targets: [^]chr|chr:pos|chr:from-to|chr:from-[,…]: Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with '^' to request logical complement. For example, '^X,Y,MT' indicates that sequences X, Y and MT should be skipped. Yet another difference between the two is that -r checks both start and end positions of indels, whereas -t checks start positions only. Note that -t cannot be used in combination with -T.
targetsFile            Optional<File>         -T                        --targets-file: Same -t, --targets, but reads regions from a file. Note that -T cannot be used in combination with -t. With the call -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele. Note that the file must be compressed and index. Such a file can be easily created from a VCF using: `bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz`
threads                Optional<Integer>      --threads                 Number of output compression threads to use in addition to main thread. Only used when --output-type is b or z. Default: 0.
siteWin                Optional<Integer>      -w                        --site-win: maximum distance between two records to consider when locally sorting variants which changed position during the realignment
=====================  =====================  ============  ==========  ============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task bcftoolsNorm {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File vcf
       String? outputFilename
       String? checkRef
       String? removeDups
       Boolean? removeDupsAcrossFiles
       File? reference
       File? reference_fai
       String? multiallelics
       Boolean? noVersion
       Boolean? noNormalize
       String? outputType
       String? regions
       File? regionsFile
       Boolean? strictFilter
       Array[File]? targets
       File? targetsFile
       Int? threads
       Int? siteWin
     }
     command <<<
       set -e
       bcftools norm \
         -o '~{select_first([outputFilename, "generated.vcf.gz"])}' \
         ~{if defined(checkRef) then ("-c '" + checkRef + "'") else ""} \
         ~{if defined(removeDups) then ("-d '" + removeDups + "'") else ""} \
         ~{if (defined(removeDupsAcrossFiles) && select_first([removeDupsAcrossFiles])) then "-D" else ""} \
         ~{if defined(reference) then ("-f '" + reference + "'") else ""} \
         ~{if defined(select_first([multiallelics, "-"])) then ("-m '" + select_first([multiallelics, "-"]) + "'") else ""} \
         ~{if (defined(noVersion) && select_first([noVersion])) then "--no-version" else ""} \
         ~{if (defined(noNormalize) && select_first([noNormalize])) then "-N" else ""} \
         ~{if defined(select_first([outputType, "z"])) then ("-O '" + select_first([outputType, "z"]) + "'") else ""} \
         ~{if defined(regions) then ("-r '" + regions + "'") else ""} \
         ~{if defined(regionsFile) then ("-R '" + regionsFile + "'") else ""} \
         ~{if (defined(strictFilter) && select_first([strictFilter])) then "-s" else ""} \
         ~{if (defined(targets) && length(select_first([targets])) > 0) then "-t '" + sep("' '", select_first([targets])) + "'" else ""} \
         ~{if defined(targetsFile) then ("-T '" + targetsFile + "'") else ""} \
         ~{if defined(threads) then ("--threads " + threads) else ''} \
         ~{if defined(siteWin) then ("-w " + siteWin) else ''} \
         '~{vcf}'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "biocontainers/bcftools:v1.5_cv2"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File out = select_first([outputFilename, "generated.vcf.gz"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: 'BCFTools: Normalize'
   doc: |
     Left-align and normalize indels, check if REF alleles match the reference,
     split multiallelic sites into multiple rows; recover multiallelics from
     multiple rows. Left-alignment and normalization will only be applied if
     the --fasta-ref option is supplied.

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: biocontainers/bcftools:v1.5_cv2

   inputs:
   - id: vcf
     label: vcf
     type: File
     inputBinding:
       position: 10
   - id: outputFilename
     label: outputFilename
     doc: |-
       --output: When output consists of a single stream, write it to FILE rather than to standard output, where it is written by default.
     type:
     - string
     - 'null'
     default: generated.vcf.gz
     inputBinding:
       prefix: -o
   - id: checkRef
     label: checkRef
     doc: |-
       --check-ref e|w|x|s: what to do when incorrect or missing REF allele is encountered: exit (e), warn (w), exclude (x), or set/fix (s) bad sites. The w option can be combined with x and s. Note that s can swap alleles and will update genotypes (GT) and AC counts, but will not attempt to fix PL or other fields. Also note, and this cannot be stressed enough, that s will NOT fix strand issues in your VCF, do NOT use it for that purpose!!! (Instead see http://samtools.github.io/bcftools/howtos/plugin.af-dist.html and http://samtools.github.io/bcftools/howtos/plugin.fixref.html.)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -c
   - id: removeDups
     label: removeDups
     doc: |-
       --rm-dup: snps|indels|both|all|none. If a record is present multiple times, output only the first instance, see --collapse in Common Options.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -d
   - id: removeDupsAcrossFiles
     label: removeDupsAcrossFiles
     doc: |-
       --remove-duplicates: If a record is present in multiple files, output only the first instance. Alias for -d none, deprecated.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -D
   - id: reference
     label: reference
     doc: |-
       --fasta-ref: reference sequence. Supplying this option will turn on left-alignment and normalization, however, see also the --do-not-normalize option below.
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .fai
     inputBinding:
       prefix: -f
   - id: multiallelics
     label: multiallelics
     doc: |-
       --multiallelics -|+[snps|indels|both|any]: split multiallelic sites into biallelic records (-) or join biallelic sites into multiallelic records (+). An optional type string can follow which controls variant types which should be split or merged together: If only SNP records should be split or merged, specify snps; if both SNPs and indels should be merged separately into two records, specify both; if SNPs and indels should be merged into a single record, specify any.
     type: string
     default: '-'
     inputBinding:
       prefix: -m
   - id: noVersion
     label: noVersion
     doc: Do not append version and command line information to the output VCF header.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no-version
   - id: noNormalize
     label: noNormalize
     doc: |-
       --do-not-normalize: the -c s option can be used to fix or set the REF allele from the reference -f. The -N option will not turn on indel normalisation as the -f option normally implies
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -N
   - id: outputType
     label: outputType
     doc: |-
       --output-type b|u|z|v: Output compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion.
     type: string
     default: z
     inputBinding:
       prefix: -O
   - id: regions
     label: regions
     doc: |-
       --regions chr|chr:pos|chr:from-to|chr:from-[,…]: Comma-separated list of regions, see also -R, --regions-file. Note that -r cannot be used in combination with -R.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: -r
   - id: regionsFile
     label: regionsFile
     doc: |-
       --regions-file: Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file are: CHROM, POS, and, optionally, POS_TO, where positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and POS_TO (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the '.bed' or '.bed.gz' suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, 'chr20' is not the same as '20'. Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -R
   - id: strictFilter
     label: strictFilter
     doc: |-
       --strict-filter: when merging (-m+), merged site is PASS only if all sites being merged PASS
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: -s
   - id: targets
     label: targets
     doc: |-
       --targets: [^]chr|chr:pos|chr:from-to|chr:from-[,…]: Similar as -r, --regions, but the next position is accessed by streaming the whole VCF/BCF rather than using the tbi/csi index. Both -r and -t options can be applied simultaneously: -r uses the index to jump to a region and -t discards positions which are not in the targets. Unlike -r, targets can be prefixed with '^' to request logical complement. For example, '^X,Y,MT' indicates that sequences X, Y and MT should be skipped. Yet another difference between the two is that -r checks both start and end positions of indels, whereas -t checks start positions only. Note that -t cannot be used in combination with -T. 
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: -t
   - id: targetsFile
     label: targetsFile
     doc: |-
       --targets-file: Same -t, --targets, but reads regions from a file. Note that -T cannot be used in combination with -t. With the call -C alleles command, third column of the targets file must be comma-separated list of alleles, starting with the reference allele. Note that the file must be compressed and index. Such a file can be easily created from a VCF using: `bcftools query -f'%CHROM\t%POS\t%REF,%ALT\n' file.vcf | bgzip -c > als.tsv.gz && tabix -s1 -b2 -e2 als.tsv.gz`
     type:
     - File
     - 'null'
     inputBinding:
       prefix: -T
   - id: threads
     label: threads
     doc: |-
       Number of output compression threads to use in addition to main thread. Only used when --output-type is b or z. Default: 0.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --threads
   - id: siteWin
     label: siteWin
     doc: |-
       --site-win: maximum distance between two records to consider when locally sorting variants which changed position during the realignment
     type:
     - int
     - 'null'
     inputBinding:
       prefix: -w

   outputs:
   - id: out
     label: out
     type: File
     outputBinding:
       glob: generated.vcf.gz
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand:
   - bcftools
   - norm
   arguments: []

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: bcftoolsNorm


