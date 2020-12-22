:orphan:

STAR Aligner
==================================

``star_genomeGenerate`` · *2 contributors · 3 versions*

Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2019 

For more details see:

- https://www.ncbi.nlm.nih.gov/pubmed/23104886
- https://github.com/alexdobin/STAR
- https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf



Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.star.versions import StarGenerateIndexes_2_5_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "star_genomegenerate_step",
           StarGenerateIndexes_2_5_3(
               outFileNamePrefix=None,
           )
       )
       wf.output("out", source=star_genomegenerate_step.out)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for star_genomeGenerate:

.. code-block:: bash

   # user inputs
   janis inputs star_genomeGenerate > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       {}




5. Run star_genomeGenerate with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       star_genomeGenerate





Information
------------

:ID: ``star_genomeGenerate``
:URL: `https://github.com/alexdobin/STAR <https://github.com/alexdobin/STAR>`_
:Versions: v2.7.5c, v2.7.1a, v2.5.3a
:Container: quay.io/biocontainers/star:2.5.3a--0
:Authors: Jiaan Yu, Michael Franklin
:Citations: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. Bioinformatics. 2013;29(1):15‐21. doi:10.1093/bioinformatics/bts635
:DOI: https://doi.org/10.1093/bioinformatics/bts635
:Created: 2020-05-29
:Updated: 2020-05-29


Outputs
-----------

======  =========  ===============
name    type       documentation
======  =========  ===============
out     Directory
======  =========  ===============


Additional configuration (inputs)
---------------------------------

================================  ========================  ==================================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                              type                      prefix                              position    documentation
================================  ========================  ==================================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
outFileNamePrefix                 String                    --outFileNamePrefix                             (default: ./) output files name prefix (including full or relative path). Can only be defined on the command line.
parametersFiles                   Optional<String>          --parametersFiles                               (default: -) none. Can only be defined on the command line.
sysShell                          Optional<String>          --sysShell                                      (default: -) path to the shell binary, preferably bash, e.g. /bin/bash.
                                                                                                            - ... the default shell is executed, typically /bin/sh. This was reported to fail on some Ubuntu systems - then you need to specify path to bash.
runThreadN                        Optional<Integer>         --runThreadN                                    (default: 1) number of threads to run STAR
runDirPerm                        Optional<String>          --runDirPerm                                    (default: User_RWX) permissions for the directories created at the run-time.
                                                                                                            - User_RWX ... user-read/write/execute
                                                                                                            - All_RWX  ... all-read/write/execute (same as chmod 777)
runRNGseed                        Optional<Integer>         --runRNGseed                                    (default: 777) random number generator seed.
genomeDir                         Optional<Directory>       --genomeDir                                     (default: GenomeDir/) path to the directory where genome files are stored
outputGenomeDir                   Optional<String>          --genomeDir                                     generated for --runMode generateGenome
genomeLoad                        Optional<String>          --genomeLoad                                    (default: NoSharedMemory) mode of shared memory usage for the genome files. Only used with --runMode alignReads.
                                                                                                            - LoadAndKeep     ... load genome into shared and keep it in memory after run,
                                                                                                            - LoadAndRemove   ... load genome into shared but remove it after run,
                                                                                                            - LoadAndExit     ... load genome into shared memory and exit, keeping the genome in memory for future runs,
                                                                                                            - Remove:   ... do not map anything, just remove loaded genome from memory,
                                                                                                            - NoSharedMemory  ... do not use shared memory, each job will have its own private copy of the genome
genomeFastaFiles                  Optional<Array<Fasta>>    --genomeFastaFiles                              (default: -) path(s) to the fasta files with the genome sequences, separated by spaces. These files should be plain text FASTA files, they *cannot* be zipped. Required for the genome generation (--runMode genomeGenerate). Can also be used in the mapping (--runMode alignReads) to add extra (new) sequences to the genome (e.g. spike-ins).
genomeChainFiles                  Optional<Array<File>>     --genomeChainFiles                              (default: -) chain files for genomic liftover. Only used with --runMode liftOver .
genomeFileSizes                   Optional<Integer>         --genomeFileSizes                               (default: 0) genome files exact sizes in bytes. Typically, this should not be defined by the user.
genomeConsensusFile               Optional<VCF>             --genomeConsensusFile                           (default: -) VCF file with consensus SNPs (i.e. alternative allele is the major (AF>0.5) allele)
genomeChrBinNbits                 Optional<Integer>         --genomeChrBinNbits                             (default: 18) each chromosome will occupy an integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as ``min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)])``.
genomeSAindexNbases               Optional<Integer>         --genomeSAindexNbases                           (default: 14) length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).
genomeSAsparseD                   Optional<Integer>         --genomeSAsparseD                               (default: 1) use bigger numbers to decrease needed RAM at the cost of mapping speed reduction
genomeSuffixLengthMax             Optional<Integer>         --genomeSuffixLengthMax                         (default: -1) maximum length of the suffixes, has to be longer than read length. -1 = infinite.
sjdbFileChrStartEnd               Optional<File>            --sjdbFileChrStartEnd                           (default: -) path to the files with genomic coordinates (chr <tab> start <tab> end <tab> strand) for the splice junction introns. Multiple files can be supplied wand will be concatenated.
sjdbGTFfile                       Optional<File>            --sjdbGTFfile                                   (default: -) path to the GTF file with annotations
sjdbGTFchrPrefix                  Optional<String>          --sjdbGTFchrPrefix                              (default: -) prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSMEBL annotations with UCSC genomes)
sjdbGTFfeatureExon                Optional<String>          --sjdbGTFfeatureExon                            (default: exon) feature type in GTF file to be used as exons for building transcripts
sjdbGTFtagExonParentTranscript    Optional<String>          --sjdbGTFtagExonParentTranscript                (default: transcript_id) GTF attribute name for parent transcript ID (default "transcript_id" works for GTF files)
sjdbGTFtagExonParentGene          Optional<String>          --sjdbGTFtagExonParentGene                      (default: gene_id) GTF attribute name for parent gene ID (default "gene_id" works for GTF files)
sjdbGTFtagExonParentGeneName      Optional<String>          --sjdbGTFtagExonParentGeneName                  (default: gene_name) GTF attrbute name for parent gene name
sjdbGTFtagExonParentGeneType      Optional<String>          --sjdbGTFtagExonParentGeneType                  (default: gene_type gene_biotype) GTF attrbute name for parent gene type
sjdbOverhang                      Optional<Integer>         --sjdbOverhang                                  (default: 100) length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
sjdbScore                         Optional<Integer>         --sjdbScore                                     (default: 2) extra alignment score for alignmets that cross database junctions
sjdbInsertSave                    Optional<String>          --sjdbInsertSave                                (default: Basic) which files to save when sjdb junctions are inserted on the fly at the mapping step
                                                                                                            - Basic ... only small junction / transcript files
                                                                                                            - All   ... all files including big Genome, SA and SAindex - this will create a complete genome directory
varVCFfile                        Optional<VCF>             --varVCFfile                                    (default: -) path to the VCF file that contains variation data.
inputBAMfile                      Optional<BAM>             --inputBAMfile                                  (default: -) path to BAM input file, to be used with --runMode inputAlignmentsFromBAM
readFilesType                     Optional<String>          --readFilesType                                 (default: Fastx) format of input read files
                                                                                                            - Fastx       ... FASTA or FASTQ
                                                                                                            - SAM SE      ... SAM or BAM single-end reads; for BAM use --readFilesCommand samtools view
                                                                                                            - SAM PE      ... SAM or BAM paired-end reads; for BAM use --readFilesCommand samtools view
readFilesIn                       Optional<FastqGzPair>     --readFilesIn                                   (default: Read1 Read2) paths to files that contain input read1 (and, if needed,  read2)
readFilesPrefix                   Optional<String>          --readFilesPrefix                               (default: -)   for the read files names, i.e. it will be added in front of the strings in --readFilesIn no prefix
readFilesCommand                  Optional<String>          --readFilesCommand                              (default: -) command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.
readMapNumber                     Optional<Integer>         --readMapNumber                                 (default: 1) number of reads to map from the beginning of the file map all reads
readMatesLengthsIn                Optional<String>          --readMatesLengthsIn                            (default: NotEqual) Equal/NotEqual - lengths of names,sequences,qualities for both mates are the same  / not the same. NotEqual is safe in all situations.
readNameSeparator                 Optional<String>          --readNameSeparator                             (default: /) character(s) separating the part of the read names that will be trimmed in output (read name after space is always trimmed)
readQualityScoreBase              Optional<Integer>         --readQualityScoreBase                          (default: 33) number to be subtracted from the ASCII code to get Phred quality score
clip3pNbases                      Optional<Integer>         --clip3pNbases                                  (default: 0) number(s) of bases to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.
clip5pNbases                      Optional<Integer>         --clip5pNbases                                  (default: 0) number(s) of bases to clip from 5p of each mate. If one value is given, it will be assumed the same for both mates.
clip3pAdapterSeq                  Optional<String>          --clip3pAdapterSeq                              (default: -) adapter sequences to clip from 3p of each mate.  If one value is given, it will be assumed the same for both mates.
clip3pAdapterMMp                  Optional<Double>          --clip3pAdapterMMp                              (default: 0.1) max proportion of mismatches for 3p adpater clipping for each mate.  If one value is given, it will be assumed the same for both mates.
clip3pAfterAdapterNbases          Optional<Integer>         --clip3pAfterAdapterNbases                      (default: 0) number of bases to clip from 3p of each mate after the adapter clipping. If one value is given, it will be assumed the same for both mates.
limitGenomeGenerateRAM            Optional<Integer>         --limitGenomeGenerateRAM                        (default: 31000000000) maximum available RAM (bytes) for genome generation
limitIObufferSize                 Optional<Integer>         --limitIObufferSize                             (default: 150000000) max available buffers size (bytes) for input/output, per thread
limitOutSAMoneReadBytes           Optional<Integer>         --limitOutSAMoneReadBytes                       (default: 100000) >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax
limitOutSJoneRead                 Optional<Integer>         --limitOutSJoneRead                             (default: 1000) max number of junctions for one read (including all multi-mappers)
limitOutSJcollapsed               Optional<Integer>         --limitOutSJcollapsed                           (default: 1000000) max number of collapsed junctions
limitBAMsortRAM                   Optional<Integer>         --limitBAMsortRAM                               (default: 0) maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with --genomeLoad NoSharedMemory option.
limitSjdbInsertNsj                Optional<Integer>         --limitSjdbInsertNsj                            (default: 1000000) maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run
limitNreadsSoft                   Optional<Integer>         --limitNreadsSoft                               (default: 1) soft limit on the number of reads
outTmpDir                         Optional<String>          --outTmpDir                                     (default: -) path to a directory that will be used as temporary by STAR. All contents of this directory will be removed!     - the temp directory will default to outFileNamePrefix_STARtmp
outTmpKeep                        Optional<String>          --outTmpKeep                                    (default: None) whether to keep the tempporary files after STAR runs is finished None ... remove all temporary files All .. keep all files
outStd                            Optional<String>          --outStd                                        (default: Log) which output will be directed to stdout (standard out) Log     ... log messages SAM        ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out BAM_Unsorted     ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate BAM_Quant        ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM
outReadsUnmapped                  Optional<String>          --outReadsUnmapped                              (default: None) which output will be directed to stdout (standard out) [Log ... log messages SAM ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out BAM_Unsorted           ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate BAM_Quant ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM]
outQSconversionAdd                Optional<Integer>         --outQSconversionAdd                            (default: 0) add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)
outMultimapperOrder               Optional<String>          --outMultimapperOrder                           (default: Old_2.4) order of multimapping alignments in the output files Old_2.4
outSAMtype                        Optional<Array<String>>   --outSAMtype                                    (default: SAM) ... quasi-random order used before 2.5.0 Random ... random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases. ... standard unsorted SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.
outSAMmode                        Optional<String>          --outSAMmode                                    (default: Full) mode of SAM output None ... no SAM output Full ... full SAM output NoQS ... full SAM but without quality scores ... no attributes Standard    ... NH HI AS nM All   ... NH HI AS nM NM MD jM jI MC ch vA    ... variant allele vG    ... genomic coordiante of the variant overlapped by the read vW    ... 0/1 - alignment does not pass / passes WASP filtering. Requires --waspOutputMode SAMtag STARsolo: CR CY UR UY ... sequences and quality scores of cell barcodes and UMIs for the solo* demultiplexing CB UB       ... error-corrected cell barcodes and UMIs for solo* demultiplexing. Requires --outSAMtype BAM SortedByCoordinate. sM    ... assessment of CB and UMI sS    ... sequence of the entire barcode (CB,UMI,adapter...) sQ    ... quality of the entire barcode Unsupported/undocumented: rB    ... alignment block read/genomic coordinates vR    ... read coordinate of the variant
outSAMstrandField                 Optional<String>          --outSAMstrandField                             (default: None) Cufflinks-like strand field flag None
outSAMattributes                  Optional<String>          --outSAMattributes                              (default: Standard) a string of desired SAM attributes, in the order desired for the output SAM NH HI AS nM NM MD jM jI XS MC ch ... any combination in any order None
outSAMattrIHstart                 Optional<Integer>         --outSAMattrIHstart                             (default: 1) start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.
outSAMunmapped                    Optional<String>          --outSAMunmapped                                (default: None) output of unmapped reads in the SAM format 1st word: None   ... no output Within ... output unmapped reads within the main SAM file (i.e. Aligned.out.sam) 2nd word: KeepPairs ... record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects multi-mapping reads.
outSAMorder                       Optional<String>          --outSAMorder                                   (default: Paired) type of sorting for the SAM output one mate after the other for all paired alignments one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files
outSAMprimaryFlag                 Optional<String>          --outSAMprimaryFlag                             (default: OneBestScore) which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG OneBestScore ... only one alignment with the best score is primary AllBestScore ... all alignments with the best score are primary
outSAMreadID                      Optional<String>          --outSAMreadID                                  (default: Standard) read ID record type Standard ... first word (until space) from the FASTx read ID line, removing /1,/2 from the end Number   ... read number (index) in the FASTx file
outSAMmapqUnique                  Optional<Integer>         --outSAMmapqUnique                              (default: 255) the MAPQ value for unique mappers
outSAMflagOR                      Optional<Integer>         --outSAMflagOR                                  (default: 0) sam FLAG will be bitwise OR'd with this value, i.e. FLAG=FLAG | outSAMflagOR. This is applied after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set otherwise.
outSAMflagAND                     Optional<Integer>         --outSAMflagAND                                 (default: 65535) sam FLAG will be bitwise AND'd with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not set otherwise.
outSAMattrRGline                  Optional<String>          --outSAMattrRGline                              (default: -) SAM/BAM read group line. The first word contains the read group identifier and must start with "ID:", e.g. --outSAMattrRGline ID:xxx CN:yy "DS:z z z".     xxx will be added as RG tag to each output alignment. Any spaces in the tag values have to be double quoted.     Comma separated RG lines correspons to different (comma separated) input files in --readFilesIn. Commas have to be surrounded by spaces, e.g.     --outSAMattrRGline ID:xxx , ID:zzz "DS:z z" , ID:yyy DS:yyyy
outSAMheaderHD                    Optional<Array<String>>   --outSAMheaderHD                                (default: -) @HD (header) line of the SAM header
outSAMheaderPG                    Optional<Array<String>>   --outSAMheaderPG                                (default: -) extra @PG (software) line of the SAM header (in addition to STAR)
outSAMheaderCommentFile           Optional<String>          --outSAMheaderCommentFile                       (default: -) path to the file with @CO (comment) lines of the SAM header
outSAMfilter                      Optional<String>          --outSAMfilter                                  (default: None) filter the output into main SAM/BAM files KeepOnlyAddedReferences ... only keep the reads for which all alignments are to the extra reference sequences added with --genomeFastaFiles at the mapping stage. KeepAllAddedReferences ...  keep all alignments to the extra reference sequences added with --genomeFastaFiles at the mapping stage.
outSAMmultNmax                    Optional<Integer>         --outSAMmultNmax                                (default: 1) max number of multiple alignments for a read that will be output to the SAM/BAM files. -1 ... all alignments (up to --outFilterMultimapNmax) will be output
outSAMtlen                        Optional<Integer>         --outSAMtlen                                    (default: 1) calculation method for the TLEN field in the SAM/BAM files 1 ... leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate 2 ... leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from 1 for overlapping mates with protruding ends
outBAMcompression                 Optional<Integer>         --outBAMcompression                             (default: 1) -1 to 10  BAM compression level, -1=default compression (6?), 0=no compression, 10=maximum compression
outBAMsortingThreadN              Optional<Integer>         --outBAMsortingThreadN                          (default: 0) number of threads for BAM sorting. 0 will default to min(6,--runThreadN).
outBAMsortingBinsN                Optional<Integer>         --outBAMsortingBinsN                            (default: 50) number of genome bins fo coordinate-sorting
bamRemoveDuplicatesType           Optional<String>          --bamRemoveDuplicatesType                       (default: -) mark duplicates in the BAM file, for now only works with (i) sorted BAM fed with inputBAMfile, and (ii) for paired-end alignments only -
bamRemoveDuplicatesMate2basesN    Optional<Integer>         --bamRemoveDuplicatesMate2basesN                (default: 0) number of bases from the 5' of mate 2 to use in collapsing (e.g. for RAMPAGE)
outWigType                        Optional<String>          --outWigType                                    (default: None) --outSAMtype BAM SortedByCoordinate .     1st word:     None       ... no signal output     bedGraph   ... bedGraph format     wiggle     ... wiggle format     2nd word:     read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc     read2      ... signal from only 2nd read
outWigStrand                      Optional<String>          --outWigStrand                                  (default: Stranded) strandedness of wiggle/bedGraph output     Stranded   ...  separate strands, str1 and str2     Unstranded ...  collapsed strands
outWigReferencesPrefix            Optional<String>          --outWigReferencesPrefix                        (default: -) prefix matching reference names to include in the output wiggle file, e.g. "chr", default "-" - include all references
outWigNorm                        Optional<String>          --outWigNorm                                    (default: RPM) type of normalization for the signal RPM    ... reads per million of mapped reads None   ... no normalization, "raw" counts
outFilterType                     Optional<String>          --outFilterType                                 (default: Normal) type of filtering Normal  ... standard filtering using only current alignment BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab
outFilterMultimapScoreRange       Optional<Integer>         --outFilterMultimapScoreRange                   (default: 1) the score range below the maximum score for multimapping alignments
outFilterMultimapNmax             Optional<Integer>         --outFilterMultimapNmax                         (default: 10) maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value.  Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out .
outFilterMismatchNmax             Optional<Integer>         --outFilterMismatchNmax                         (default: 10) alignment will be output only if it has no more mismatches than this value.
outFilterMismatchNoverLmax        Optional<Float>           --outFilterMismatchNoverLmax                    (default: 0.3) alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.
outFilterMismatchNoverReadLmax    Optional<Float>           --outFilterMismatchNoverReadLmax                (default: 1) alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value.
outFilterScoreMin                 Optional<Integer>         --outFilterScoreMin                             (default: 0) alignment will be output only if its score is higher than or equal to this value.
outFilterScoreMinOverLread        Optional<Float>           --outFilterScoreMinOverLread                    (default: 0.66) same as outFilterScoreMin, but  normalized to read length (sum of mates' lengths for paired-end reads)
outFilterMatchNmin                Optional<Integer>         --outFilterMatchNmin                            (default: 0) alignment will be output only if the number of matched bases is higher than or equal to this value.
outFilterMatchNminOverLread       Optional<Float>           --outFilterMatchNminOverLread                   (default: 0.66) sam as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads).
outFilterIntronMotifs             Optional<String>          --outFilterIntronMotifs                         (default: None) filter alignment using their motifs None
outFilterIntronStrands            Optional<String>          --outFilterIntronStrands                        (default: RemoveInconsistentStrands) filter alignments RemoveInconsistentStrands      ... remove alignments that have junctions with inconsistent strands None
outSJfilterReads                  Optional<String>          --outSJfilterReads                              (default: All) which reads to consider for collapsed splice junctions output all reads, unique- and multi-mappers uniquely mapping reads only
outSJfilterOverhangMin            Optional<Integer>         --outSJfilterOverhangMin                        (default: 30 12 12 12) minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif does not apply to annotated junctions
outSJfilterCountUniqueMin         Optional<Integer>         --outSJfilterCountUniqueMin                     (default: 3 1 1 1) minimum uniquely mapping read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied does not apply to annotated junctions
outSJfilterCountTotalMin          Optional<Integer>         --outSJfilterCountTotalMin                      (default: 3 1 1 1) minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied does not apply to annotated junctions
outSJfilterDistToOtherSJmin       Optional<Integer>         --outSJfilterDistToOtherSJmin                   (default: 10 0 5 10) minimum allowed distance to other junctions' donor/acceptor does not apply to annotated junctions
outSJfilterIntronMaxVsReadN       Optional<Integer>         --outSJfilterIntronMaxVsReadN                   (default: 50000 100000 200000) maximum gap allowed for junctions supported by 1,2,3,,,N reads <=200000. by >=4 reads any gap <=alignIntronMax does not apply to annotated junctions
scoreGap                          Optional<Integer>         --scoreGap                                      (default: 0) splice junction penalty (independent on intron motif)
scoreGapNoncan                    Optional<Integer>         --scoreGapNoncan                                (default: 8) non-canonical junction penalty (in addition to scoreGap)
scoreGapGCAG                      Optional<Integer>         --scoreGapGCAG                                  (default: 4) GC/AG and CT/GC junction penalty (in addition to scoreGap)
scoreGapATAC                      Optional<Integer>         --scoreGapATAC                                  (default: 8) AT/AC  and GT/AT junction penalty  (in addition to scoreGap)
scoreGenomicLengthLog2scale       Optional<Float>           --scoreGenomicLengthLog2scale                   (default: -0.25) scoreGenomicLengthLog2scale*log2(genomicLength)
scoreDelOpen                      Optional<Integer>         --scoreDelOpen                                  (default: 2) deletion open penalty
scoreDelBase                      Optional<Integer>         --scoreDelBase                                  (default: 2) deletion extension penalty per base (in addition to scoreDelOpen)
scoreInsOpen                      Optional<Integer>         --scoreInsOpen                                  (default: 2) insertion open penalty
scoreInsBase                      Optional<Integer>         --scoreInsBase                                  (default: 2) insertion extension penalty per base (in addition to scoreInsOpen)
scoreStitchSJshift                Optional<Integer>         --scoreStitchSJshift                            (default: 1) maximum score reduction while searching for SJ boundaries inthe stitching step
seedSearchStartLmax               Optional<Integer>         --seedSearchStartLmax                           (default: 50) defines the search start point through the read - the read is split into pieces no longer than this value
seedSearchStartLmaxOverLread      Optional<Float>           --seedSearchStartLmaxOverLread                  (default: 1) seedSearchStartLmax normalized to read length (sum of mates' lengths for paired-end reads)
seedSearchLmax                    Optional<Integer>         --seedSearchLmax                                (default: 0) defines the maximum length of the seeds, if =0 max seed lengthis infinite
seedMultimapNmax                  Optional<Integer>         --seedMultimapNmax                              (default: 10000) only pieces that map fewer than this value are utilized in the stitching procedure
seedPerReadNmax                   Optional<Integer>         --seedPerReadNmax                               (default: 1000) max number of seeds per read
seedPerWindowNmax                 Optional<Integer>         --seedPerWindowNmax                             (default: 50) max number of seeds per window
seedNoneLociPerWindow             Optional<Integer>         --seedNoneLociPerWindow                         (default: 10) max number of one seed loci per window
seedSplitMin                      Optional<Integer>         --seedSplitMin                                  (default: 12) min length of the seed sequences split by Ns or mate gap
alignIntronMin                    Optional<Integer>         --alignIntronMin                                (default: 21) genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
alignIntronMax                    Optional<Integer>         --alignIntronMax                                (default: 0) maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins
alignMatesGapMax                  Optional<Integer>         --alignMatesGapMax                              (default: 0) maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins
alignSJoverhangMin                Optional<Integer>         --alignSJoverhangMin                            (default: 5) minimum overhang (i.e. block size) for spliced alignments
alignSJstitchMismatchNmax         Optional<Array<Integer>>  --alignSJstitchMismatchNmax                     (default: 0 -1 0 0) maximum number of mismatches for stitching of the splice junctions (-1: no limit).     (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.
alignSJDBoverhangMin              Optional<Integer>         --alignSJDBoverhangMin                          (default: 3) minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments
alignSplicedMateMapLmin           Optional<Integer>         --alignSplicedMateMapLmin                       (default: 0) minimum mapped length for a read mate that is spliced
alignSplicedMateMapLminOverLmate  Optional<Float>           --alignSplicedMateMapLminOverLmate              (default: 0.66) alignSplicedMateMapLmin normalized to mate length
alignWindowsPerReadNmax           Optional<Integer>         --alignWindowsPerReadNmax                       (default: 10000) max number of windows per read
alignTranscriptsPerWindowNmax     Optional<Integer>         --alignTranscriptsPerWindowNmax                 (default: 100) max number of transcripts per window
alignTranscriptsPerReadNmax       Optional<Integer>         --alignTranscriptsPerReadNmax                   (default: 10000) max number of different alignments per read to consider
alignEndsType                     Optional<String>          --alignEndsType                                 (default: Local) type of read ends alignment Local
alignEndsProtrude                 Optional<Integer>         --alignEndsProtrude                             (default: 0 ConcordantPair) allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate maximum number of protrusion bases allowed string:     ConcordantPair ... report alignments with non-zero protrusion as concordant pairs     DiscordantPair ... report alignments with non-zero protrusion as discordant pairs
alignSoftClipAtReferenceEnds      Optional<String>          --alignSoftClipAtReferenceEnds                  (default: Yes) allow the soft-clipping of the alignments past the end of the chromosomes Yes ... allow No  ... prohibit, useful for compatibility with Cufflinks
alignInsertionFlush               Optional<String>          --alignInsertionFlush                           (default: None) how to flush ambiguous insertion positions None    ... insertions are not flushed Right   ... insertions are flushed to the right
peOverlapNbasesMin                Optional<Integer>         --peOverlapNbasesMin                            (default: 0) minimum number of overlap bases to trigger mates merging and realignment
peOverlapMMp                      Optional<Float>           --peOverlapMMp                                  (default: 0.01) maximum proportion of mismatched bases in the overlap area
winAnchorMultimapNmax             Optional<Integer>         --winAnchorMultimapNmax                         (default: 50) max number of loci anchors are allowed to map to
winBinNbits                       Optional<Integer>         --winBinNbits                                   (default: 16) =LOG2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins.
winAnchorDistNbins                Optional<Integer>         --winAnchorDistNbins                            (default: 9) max number of bins between two anchors that allows aggregation of anchors into one window
winFlankNbins                     Optional<Integer>         --winFlankNbins                                 (default: 4) log2(winFlank), where win Flank is the size of the left and right flanking regions for each window
winReadCoverageRelativeMin        Optional<Float>           --winReadCoverageRelativeMin                    (default: 0.5) minimum relative coverage of the read sequence by the seeds in a window, for STARlong algorithm only.
winReadCoverageBasesMin           Optional<Integer>         --winReadCoverageBasesMin                       (default: 0) minimum number of bases covered by the seeds in a window , for STARlong algorithm only.
chimOutType                       Optional<Array<String>>   --chimOutType                                   (default: Junctions) type of chimeric output     Junctions       ... Chimeric.out.junction     SeparateSAMold  ... output old SAM into separate Chimeric.out.sam file     WithinBAM       ... output into main aligned BAM files (Aligned.*.bam)     WithinBAM HardClip  ... (default) hard-clipping in the CIGAR for supplemental chimeric alignments (defaultif no 2nd word is present)     WithinBAM SoftClip  ... soft-clipping in the CIGAR for supplemental chimeric alignments
chimSegmentMin                    Optional<Integer>         --chimSegmentMin                                (default: 0) minimum length of chimeric segment length, if ==0, no chimeric output
chimScoreMin                      Optional<Integer>         --chimScoreMin                                  (default: 0) minimum total (summed) score of the chimeric segments
chimScoreDropMax                  Optional<Integer>         --chimScoreDropMax                              (default: 20) max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length
chimScoreSeparation               Optional<Integer>         --chimScoreSeparation                           (default: 10) minimum difference (separation) between the best chimeric score and the next one
chimScoreJunctionNonGTAG          Optional<Integer>         --chimScoreJunctionNonGTAG                      (default: -1) penalty for a non-GT/AG chimeric junction
chimJunctionOverhangMin           Optional<Integer>         --chimJunctionOverhangMin                       (default: 20) minimum overhang for a chimeric junction
chimSegmentReadGapMax             Optional<Integer>         --chimSegmentReadGapMax                         (default: 0) maximum gap in the read sequence between chimeric segments
chimFilter                        Optional<String>          --chimFilter                                    (default: banGenomicN) different filters for chimeric alignments     None ... no filtering     banGenomicN ... Ns are not allowed in the genome sequence around the chimeric junction
chimMainSegmentMultNmax           Optional<Integer>         --chimMainSegmentMultNmax                       (default: 10) maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments.
chimMultimapNmax                  Optional<Integer>         --chimMultimapNmax                              (default: 0) maximum number of chimeric multi-alignments 0 ... use the old scheme for chimeric detection which only considered unique alignments
chimMultimapScoreRange            Optional<Integer>         --chimMultimapScoreRange                        (default: 1) the score range for multi-mapping chimeras below the best chimeric score. Only works with --chimMultimapNmax > 1
chimNonchimScoreDropMin           Optional<Integer>         --chimNonchimScoreDropMin                       (default: 20) to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value ... none     TranscriptomeSAM ... output SAM/BAM alignments to transcriptome into a separate file     GeneCounts       ... count reads per gene
chimOutJunctionFormat             Optional<Integer>         --chimOutJunctionFormat                         (default: 0) formatting type for the Chimeric.out.junction file 0 ... no comment lines/headers total, unique, multi
quantMode                         Optional<String>          --quantMode                                     (default: -) types of quantification requested     -        ... prohibit single-end alignments
quantTranscriptomeBAMcompression  Optional<Integer>         --quantTranscriptomeBAMcompression              (default: 1 1) -2 to 10  transcriptome BAM compression level     -2  ... no BAM output     -1  ... default compression (6?)      0  ... no compression      10 ... maximum compression ... 1-pass mapping     Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
quantTranscriptomeBan             Optional<String>          --quantTranscriptomeBan                         (default: IndelSoftclipSingleend) prohibit various alignment type     IndelSoftclipSingleend  ... prohibit indels, soft clipping and single-end alignments - compatible with RSEM     Singleend
twopassMode                       Optional<String>          --twopassMode                                   (default: None) 2-pass mapping mode.     None
twopass1readsN                    Optional<Integer>         --twopass1readsN                                (default: 1) number of reads to process for the 1st step. Use very large number (or default -1) to map all reads in the first step.
waspOutputMode                    Optional<String>          --waspOutputMode                                (default: None) Nature Methods 12, 1061–1063 (2015), https://www.nature.com/articles/nmeth.3582 .     SAMtag      ... add WASP tags to the alignments that pass WASP filtering
soloType                          Optional<String>          --soloType                                      (default: None) type of single-cell RNA-seq     CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium     CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.
soloCBwhitelist                   Optional<String>          --soloCBwhitelist                               (default: -) file(s) with whitelist(s) of cell barcodes. Only one file allowed with
soloCBstart                       Optional<Integer>         --soloCBstart                                   (default: 1) cell barcode start base
soloCBlen                         Optional<Integer>         --soloCBlen                                     (default: 16) cell barcode length
soloUMIstart                      Optional<Integer>         --soloUMIstart                                  (default: 17) UMI start base
soloUMIlen                        Optional<Integer>         --soloUMIlen                                    (default: 10) UMI length
soloBarcodeReadLength             Optional<Integer>         --soloBarcodeReadLength                         (default: 1) length of the barcode read     1   ... equal to sum of soloCBlen+soloUMIlen     0   ... not defined, do not check
soloCBposition                    Optional<Array<String>>   --soloCBposition                                (default: -) position of Cell Barcode(s) on the barcode read.     Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2. startAnchor_startDistance_endAnchor_endDistance adapter end     start(end)Distance is the distance from the CB start(end) to the Anchor base     String for different barcodes are separated by space. inDrop (Zilionis et al, Nat. Protocols, 2017):     --soloCBposition  0_0_2_-1  3_1_3_8
soloUMIposition                   Optional<String>          --soloUMIposition                               (default: -) position of the UMI on the barcode read, same as soloCBposition inDrop (Zilionis et al, Nat. Protocols, 2017):     --soloCBposition  3_9_3_14
soloAdapterSequence               Optional<String>          --soloAdapterSequence                           (default: -) adapter sequence to anchor barcodes.    ... only exact matches allowed     1MM         ... only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.     1MM_multi         ... multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches.  Allowed CBs have to have at least one read with exact match. Similar to CellRanger 2.2.0     1MM_multi_pseudocounts  ... same as 1MM_Multi, but pseudocounts of 1 are added to all whitelist barcodes. Similar to CellRanger 3.x.x
soloAdapterMismatchesNmax         Optional<Integer>         --soloAdapterMismatchesNmax                     (default: 1) maximum number of mismatches allowed in adapter sequence.
soloCBmatchWLtype                 Optional<String>          --soloCBmatchWLtype                             (default: 1MM_multi) matching the Cell Barcodes to the WhiteList     Exact
soloStrand                        Optional<String>          --soloStrand                                    (default: Forward) strandedness of the solo libraries:     Unstranded  ... no strand information     Forward     ... read strand same as the original RNA molecule     Reverse     ... read strand opposite to the original RNA molecule .. all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)     1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).     Exact       ... only exactly matching UMIs are collapsed
soloFeatures                      Optional<String>          --soloFeatures                                  (default: Gene) genomic features for which the UMI counts per Cell Barcode are collected reads match the gene transcript reported in SJ.out.tab count all reads overlapping genes' exons and introns     Transcript3p   ... quantification of transcript for 3' protocols
soloUMIdedup                      Optional<String>          --soloUMIdedup                                  (default: 1MM_All) type of UMI deduplication (collapsing) algorithm     1MM_All
soloUMIfiltering                  Optional<String>          --soloUMIfiltering                              (default: -) type of UMI filtering remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)     MultiGeneUMI    ... remove lower-count UMIs that map to more than one gene (introduced in CellRanger 3.x.x)
soloOutFileNames                  Optional<String>          --soloOutFileNames                              (default: Solo.out/  features.tsv barcodes.tsv matrix.mtx) file names for STARsolo output:     file_name_prefix   gene_names   barcode_sequences   cell_feature_count_matrix
soloCellFilter                    Optional<String>          --soloCellFilter                                (default: CellRanger2.2 3000 0.99 10) ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)     1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).     Exact       ... only exactly matching UMIs are collapsed
================================  ========================  ==================================  ==========  ==============================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task star_genomeGenerate {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       String? parametersFiles
       String? sysShell
       Int? runThreadN
       String? runDirPerm
       Int? runRNGseed
       Directory? genomeDir
       String? outputGenomeDir
       String? genomeLoad
       Array[File]? genomeFastaFiles
       Array[File]? genomeChainFiles
       Int? genomeFileSizes
       File? genomeConsensusFile
       Int? genomeChrBinNbits
       Int? genomeSAindexNbases
       Int? genomeSAsparseD
       Int? genomeSuffixLengthMax
       File? sjdbFileChrStartEnd
       File? sjdbGTFfile
       String? sjdbGTFchrPrefix
       String? sjdbGTFfeatureExon
       String? sjdbGTFtagExonParentTranscript
       String? sjdbGTFtagExonParentGene
       String? sjdbGTFtagExonParentGeneName
       String? sjdbGTFtagExonParentGeneType
       Int? sjdbOverhang
       Int? sjdbScore
       String? sjdbInsertSave
       File? varVCFfile
       File? inputBAMfile
       String? readFilesType
       Array[File]? readFilesIn
       String? readFilesPrefix
       String? readFilesCommand
       Int? readMapNumber
       String? readMatesLengthsIn
       String? readNameSeparator
       Int? readQualityScoreBase
       Int? clip3pNbases
       Int? clip5pNbases
       String? clip3pAdapterSeq
       Float? clip3pAdapterMMp
       Int? clip3pAfterAdapterNbases
       Int? limitGenomeGenerateRAM
       Int? limitIObufferSize
       Int? limitOutSAMoneReadBytes
       Int? limitOutSJoneRead
       Int? limitOutSJcollapsed
       Int? limitBAMsortRAM
       Int? limitSjdbInsertNsj
       Int? limitNreadsSoft
       String? outFileNamePrefix
       String? outTmpDir
       String? outTmpKeep
       String? outStd
       String? outReadsUnmapped
       Int? outQSconversionAdd
       String? outMultimapperOrder
       Array[String]? outSAMtype
       String? outSAMmode
       String? outSAMstrandField
       String? outSAMattributes
       Int? outSAMattrIHstart
       String? outSAMunmapped
       String? outSAMorder
       String? outSAMprimaryFlag
       String? outSAMreadID
       Int? outSAMmapqUnique
       Int? outSAMflagOR
       Int? outSAMflagAND
       String? outSAMattrRGline
       Array[String]? outSAMheaderHD
       Array[String]? outSAMheaderPG
       String? outSAMheaderCommentFile
       String? outSAMfilter
       Int? outSAMmultNmax
       Int? outSAMtlen
       Int? outBAMcompression
       Int? outBAMsortingThreadN
       Int? outBAMsortingBinsN
       String? bamRemoveDuplicatesType
       Int? bamRemoveDuplicatesMate2basesN
       String? outWigType
       String? outWigStrand
       String? outWigReferencesPrefix
       String? outWigNorm
       String? outFilterType
       Int? outFilterMultimapScoreRange
       Int? outFilterMultimapNmax
       Int? outFilterMismatchNmax
       Float? outFilterMismatchNoverLmax
       Float? outFilterMismatchNoverReadLmax
       Int? outFilterScoreMin
       Float? outFilterScoreMinOverLread
       Int? outFilterMatchNmin
       Float? outFilterMatchNminOverLread
       String? outFilterIntronMotifs
       String? outFilterIntronStrands
       String? outSJfilterReads
       Int? outSJfilterOverhangMin
       Int? outSJfilterCountUniqueMin
       Int? outSJfilterCountTotalMin
       Int? outSJfilterDistToOtherSJmin
       Int? outSJfilterIntronMaxVsReadN
       Int? scoreGap
       Int? scoreGapNoncan
       Int? scoreGapGCAG
       Int? scoreGapATAC
       Float? scoreGenomicLengthLog2scale
       Int? scoreDelOpen
       Int? scoreDelBase
       Int? scoreInsOpen
       Int? scoreInsBase
       Int? scoreStitchSJshift
       Int? seedSearchStartLmax
       Float? seedSearchStartLmaxOverLread
       Int? seedSearchLmax
       Int? seedMultimapNmax
       Int? seedPerReadNmax
       Int? seedPerWindowNmax
       Int? seedNoneLociPerWindow
       Int? seedSplitMin
       Int? alignIntronMin
       Int? alignIntronMax
       Int? alignMatesGapMax
       Int? alignSJoverhangMin
       Array[Int]? alignSJstitchMismatchNmax
       Int? alignSJDBoverhangMin
       Int? alignSplicedMateMapLmin
       Float? alignSplicedMateMapLminOverLmate
       Int? alignWindowsPerReadNmax
       Int? alignTranscriptsPerWindowNmax
       Int? alignTranscriptsPerReadNmax
       String? alignEndsType
       Int? alignEndsProtrude
       String? alignSoftClipAtReferenceEnds
       String? alignInsertionFlush
       Int? peOverlapNbasesMin
       Float? peOverlapMMp
       Int? winAnchorMultimapNmax
       Int? winBinNbits
       Int? winAnchorDistNbins
       Int? winFlankNbins
       Float? winReadCoverageRelativeMin
       Int? winReadCoverageBasesMin
       Array[String]? chimOutType
       Int? chimSegmentMin
       Int? chimScoreMin
       Int? chimScoreDropMax
       Int? chimScoreSeparation
       Int? chimScoreJunctionNonGTAG
       Int? chimJunctionOverhangMin
       Int? chimSegmentReadGapMax
       String? chimFilter
       Int? chimMainSegmentMultNmax
       Int? chimMultimapNmax
       Int? chimMultimapScoreRange
       Int? chimNonchimScoreDropMin
       Int? chimOutJunctionFormat
       String? quantMode
       Int? quantTranscriptomeBAMcompression
       String? quantTranscriptomeBan
       String? twopassMode
       Int? twopass1readsN
       String? waspOutputMode
       String? soloType
       String? soloCBwhitelist
       Int? soloCBstart
       Int? soloCBlen
       Int? soloUMIstart
       Int? soloUMIlen
       Int? soloBarcodeReadLength
       Array[String]? soloCBposition
       String? soloUMIposition
       String? soloAdapterSequence
       Int? soloAdapterMismatchesNmax
       String? soloCBmatchWLtype
       String? soloStrand
       String? soloFeatures
       String? soloUMIdedup
       String? soloUMIfiltering
       String? soloOutFileNames
       String? soloCellFilter
     }
     command <<<
       set -e
       STAR \
         ~{if defined(parametersFiles) then ("--parametersFiles '" + parametersFiles + "'") else ""} \
         ~{if defined(sysShell) then ("--sysShell '" + sysShell + "'") else ""} \
         ~{if defined(select_first([runThreadN, select_first([runtime_cpu, 1])])) then ("--runThreadN " + select_first([runThreadN, select_first([runtime_cpu, 1])])) else ''} \
         ~{if defined(runDirPerm) then ("--runDirPerm '" + runDirPerm + "'") else ""} \
         ~{if defined(runRNGseed) then ("--runRNGseed " + runRNGseed) else ''} \
         ~{if defined(genomeDir) then ("--genomeDir '" + genomeDir + "'") else ""} \
         ~{if defined(outputGenomeDir) then ("--genomeDir '" + outputGenomeDir + "'") else ""} \
         ~{if defined(genomeLoad) then ("--genomeLoad '" + genomeLoad + "'") else ""} \
         ~{if (defined(genomeFastaFiles) && length(select_first([genomeFastaFiles])) > 0) then "--genomeFastaFiles '" + sep("' '", select_first([genomeFastaFiles])) + "'" else ""} \
         ~{if (defined(genomeChainFiles) && length(select_first([genomeChainFiles])) > 0) then "--genomeChainFiles '" + sep("' '", select_first([genomeChainFiles])) + "'" else ""} \
         ~{if defined(genomeFileSizes) then ("--genomeFileSizes " + genomeFileSizes) else ''} \
         ~{if defined(genomeConsensusFile) then ("--genomeConsensusFile '" + genomeConsensusFile + "'") else ""} \
         ~{if defined(genomeChrBinNbits) then ("--genomeChrBinNbits " + genomeChrBinNbits) else ''} \
         ~{if defined(genomeSAindexNbases) then ("--genomeSAindexNbases " + genomeSAindexNbases) else ''} \
         ~{if defined(genomeSAsparseD) then ("--genomeSAsparseD " + genomeSAsparseD) else ''} \
         ~{if defined(genomeSuffixLengthMax) then ("--genomeSuffixLengthMax " + genomeSuffixLengthMax) else ''} \
         ~{if defined(sjdbFileChrStartEnd) then ("--sjdbFileChrStartEnd '" + sjdbFileChrStartEnd + "'") else ""} \
         ~{if defined(sjdbGTFfile) then ("--sjdbGTFfile '" + sjdbGTFfile + "'") else ""} \
         ~{if defined(sjdbGTFchrPrefix) then ("--sjdbGTFchrPrefix '" + sjdbGTFchrPrefix + "'") else ""} \
         ~{if defined(sjdbGTFfeatureExon) then ("--sjdbGTFfeatureExon '" + sjdbGTFfeatureExon + "'") else ""} \
         ~{if defined(sjdbGTFtagExonParentTranscript) then ("--sjdbGTFtagExonParentTranscript '" + sjdbGTFtagExonParentTranscript + "'") else ""} \
         ~{if defined(sjdbGTFtagExonParentGene) then ("--sjdbGTFtagExonParentGene '" + sjdbGTFtagExonParentGene + "'") else ""} \
         ~{if defined(sjdbGTFtagExonParentGeneName) then ("--sjdbGTFtagExonParentGeneName '" + sjdbGTFtagExonParentGeneName + "'") else ""} \
         ~{if defined(sjdbGTFtagExonParentGeneType) then ("--sjdbGTFtagExonParentGeneType '" + sjdbGTFtagExonParentGeneType + "'") else ""} \
         ~{if defined(sjdbOverhang) then ("--sjdbOverhang " + sjdbOverhang) else ''} \
         ~{if defined(sjdbScore) then ("--sjdbScore " + sjdbScore) else ''} \
         ~{if defined(sjdbInsertSave) then ("--sjdbInsertSave '" + sjdbInsertSave + "'") else ""} \
         ~{if defined(varVCFfile) then ("--varVCFfile '" + varVCFfile + "'") else ""} \
         ~{if defined(inputBAMfile) then ("--inputBAMfile '" + inputBAMfile + "'") else ""} \
         ~{if defined(readFilesType) then ("--readFilesType '" + readFilesType + "'") else ""} \
         ~{if (defined(readFilesIn) && length(select_first([readFilesIn])) > 0) then "--readFilesIn '" + sep("' '", select_first([readFilesIn])) + "'" else ""} \
         ~{if defined(readFilesPrefix) then ("--readFilesPrefix '" + readFilesPrefix + "'") else ""} \
         ~{if defined(readFilesCommand) then ("--readFilesCommand '" + readFilesCommand + "'") else ""} \
         ~{if defined(readMapNumber) then ("--readMapNumber " + readMapNumber) else ''} \
         ~{if defined(readMatesLengthsIn) then ("--readMatesLengthsIn '" + readMatesLengthsIn + "'") else ""} \
         ~{if defined(readNameSeparator) then ("--readNameSeparator '" + readNameSeparator + "'") else ""} \
         ~{if defined(readQualityScoreBase) then ("--readQualityScoreBase " + readQualityScoreBase) else ''} \
         ~{if defined(clip3pNbases) then ("--clip3pNbases " + clip3pNbases) else ''} \
         ~{if defined(clip5pNbases) then ("--clip5pNbases " + clip5pNbases) else ''} \
         ~{if defined(clip3pAdapterSeq) then ("--clip3pAdapterSeq '" + clip3pAdapterSeq + "'") else ""} \
         ~{if defined(clip3pAdapterMMp) then ("--clip3pAdapterMMp " + clip3pAdapterMMp) else ''} \
         ~{if defined(clip3pAfterAdapterNbases) then ("--clip3pAfterAdapterNbases " + clip3pAfterAdapterNbases) else ''} \
         ~{if defined(limitGenomeGenerateRAM) then ("--limitGenomeGenerateRAM " + limitGenomeGenerateRAM) else ''} \
         ~{if defined(limitIObufferSize) then ("--limitIObufferSize " + limitIObufferSize) else ''} \
         ~{if defined(limitOutSAMoneReadBytes) then ("--limitOutSAMoneReadBytes " + limitOutSAMoneReadBytes) else ''} \
         ~{if defined(limitOutSJoneRead) then ("--limitOutSJoneRead " + limitOutSJoneRead) else ''} \
         ~{if defined(limitOutSJcollapsed) then ("--limitOutSJcollapsed " + limitOutSJcollapsed) else ''} \
         ~{if defined(limitBAMsortRAM) then ("--limitBAMsortRAM " + limitBAMsortRAM) else ''} \
         ~{if defined(limitSjdbInsertNsj) then ("--limitSjdbInsertNsj " + limitSjdbInsertNsj) else ''} \
         ~{if defined(limitNreadsSoft) then ("--limitNreadsSoft " + limitNreadsSoft) else ''} \
         --outFileNamePrefix '~{select_first([outFileNamePrefix, "./"])}' \
         ~{if defined(outTmpDir) then ("--outTmpDir '" + outTmpDir + "'") else ""} \
         ~{if defined(outTmpKeep) then ("--outTmpKeep '" + outTmpKeep + "'") else ""} \
         ~{if defined(outStd) then ("--outStd '" + outStd + "'") else ""} \
         ~{if defined(outReadsUnmapped) then ("--outReadsUnmapped '" + outReadsUnmapped + "'") else ""} \
         ~{if defined(outQSconversionAdd) then ("--outQSconversionAdd " + outQSconversionAdd) else ''} \
         ~{if defined(outMultimapperOrder) then ("--outMultimapperOrder '" + outMultimapperOrder + "'") else ""} \
         ~{if (defined(outSAMtype) && length(select_first([outSAMtype])) > 0) then "--outSAMtype '" + sep("' '", select_first([outSAMtype])) + "'" else ""} \
         ~{if defined(outSAMmode) then ("--outSAMmode '" + outSAMmode + "'") else ""} \
         ~{if defined(outSAMstrandField) then ("--outSAMstrandField '" + outSAMstrandField + "'") else ""} \
         ~{if defined(outSAMattributes) then ("--outSAMattributes '" + outSAMattributes + "'") else ""} \
         ~{if defined(outSAMattrIHstart) then ("--outSAMattrIHstart " + outSAMattrIHstart) else ''} \
         ~{if defined(outSAMunmapped) then ("--outSAMunmapped '" + outSAMunmapped + "'") else ""} \
         ~{if defined(outSAMorder) then ("--outSAMorder '" + outSAMorder + "'") else ""} \
         ~{if defined(outSAMprimaryFlag) then ("--outSAMprimaryFlag '" + outSAMprimaryFlag + "'") else ""} \
         ~{if defined(outSAMreadID) then ("--outSAMreadID '" + outSAMreadID + "'") else ""} \
         ~{if defined(outSAMmapqUnique) then ("--outSAMmapqUnique " + outSAMmapqUnique) else ''} \
         ~{if defined(outSAMflagOR) then ("--outSAMflagOR " + outSAMflagOR) else ''} \
         ~{if defined(outSAMflagAND) then ("--outSAMflagAND " + outSAMflagAND) else ''} \
         ~{if defined(outSAMattrRGline) then ("--outSAMattrRGline '" + outSAMattrRGline + "'") else ""} \
         ~{if (defined(outSAMheaderHD) && length(select_first([outSAMheaderHD])) > 0) then "--outSAMheaderHD '" + sep("' '", select_first([outSAMheaderHD])) + "'" else ""} \
         ~{if (defined(outSAMheaderPG) && length(select_first([outSAMheaderPG])) > 0) then "--outSAMheaderPG '" + sep("' '", select_first([outSAMheaderPG])) + "'" else ""} \
         ~{if defined(outSAMheaderCommentFile) then ("--outSAMheaderCommentFile '" + outSAMheaderCommentFile + "'") else ""} \
         ~{if defined(outSAMfilter) then ("--outSAMfilter '" + outSAMfilter + "'") else ""} \
         ~{if defined(outSAMmultNmax) then ("--outSAMmultNmax " + outSAMmultNmax) else ''} \
         ~{if defined(outSAMtlen) then ("--outSAMtlen " + outSAMtlen) else ''} \
         ~{if defined(outBAMcompression) then ("--outBAMcompression " + outBAMcompression) else ''} \
         ~{if defined(outBAMsortingThreadN) then ("--outBAMsortingThreadN " + outBAMsortingThreadN) else ''} \
         ~{if defined(outBAMsortingBinsN) then ("--outBAMsortingBinsN " + outBAMsortingBinsN) else ''} \
         ~{if defined(bamRemoveDuplicatesType) then ("--bamRemoveDuplicatesType '" + bamRemoveDuplicatesType + "'") else ""} \
         ~{if defined(bamRemoveDuplicatesMate2basesN) then ("--bamRemoveDuplicatesMate2basesN " + bamRemoveDuplicatesMate2basesN) else ''} \
         ~{if defined(outWigType) then ("--outWigType '" + outWigType + "'") else ""} \
         ~{if defined(outWigStrand) then ("--outWigStrand '" + outWigStrand + "'") else ""} \
         ~{if defined(outWigReferencesPrefix) then ("--outWigReferencesPrefix '" + outWigReferencesPrefix + "'") else ""} \
         ~{if defined(outWigNorm) then ("--outWigNorm '" + outWigNorm + "'") else ""} \
         ~{if defined(outFilterType) then ("--outFilterType '" + outFilterType + "'") else ""} \
         ~{if defined(outFilterMultimapScoreRange) then ("--outFilterMultimapScoreRange " + outFilterMultimapScoreRange) else ''} \
         ~{if defined(outFilterMultimapNmax) then ("--outFilterMultimapNmax " + outFilterMultimapNmax) else ''} \
         ~{if defined(outFilterMismatchNmax) then ("--outFilterMismatchNmax " + outFilterMismatchNmax) else ''} \
         ~{if defined(outFilterMismatchNoverLmax) then ("--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax) else ''} \
         ~{if defined(outFilterMismatchNoverReadLmax) then ("--outFilterMismatchNoverReadLmax " + outFilterMismatchNoverReadLmax) else ''} \
         ~{if defined(outFilterScoreMin) then ("--outFilterScoreMin " + outFilterScoreMin) else ''} \
         ~{if defined(outFilterScoreMinOverLread) then ("--outFilterScoreMinOverLread " + outFilterScoreMinOverLread) else ''} \
         ~{if defined(outFilterMatchNmin) then ("--outFilterMatchNmin " + outFilterMatchNmin) else ''} \
         ~{if defined(outFilterMatchNminOverLread) then ("--outFilterMatchNminOverLread " + outFilterMatchNminOverLread) else ''} \
         ~{if defined(outFilterIntronMotifs) then ("--outFilterIntronMotifs '" + outFilterIntronMotifs + "'") else ""} \
         ~{if defined(outFilterIntronStrands) then ("--outFilterIntronStrands '" + outFilterIntronStrands + "'") else ""} \
         ~{if defined(outSJfilterReads) then ("--outSJfilterReads '" + outSJfilterReads + "'") else ""} \
         ~{if defined(outSJfilterOverhangMin) then ("--outSJfilterOverhangMin " + outSJfilterOverhangMin) else ''} \
         ~{if defined(outSJfilterCountUniqueMin) then ("--outSJfilterCountUniqueMin " + outSJfilterCountUniqueMin) else ''} \
         ~{if defined(outSJfilterCountTotalMin) then ("--outSJfilterCountTotalMin " + outSJfilterCountTotalMin) else ''} \
         ~{if defined(outSJfilterDistToOtherSJmin) then ("--outSJfilterDistToOtherSJmin " + outSJfilterDistToOtherSJmin) else ''} \
         ~{if defined(outSJfilterIntronMaxVsReadN) then ("--outSJfilterIntronMaxVsReadN " + outSJfilterIntronMaxVsReadN) else ''} \
         ~{if defined(scoreGap) then ("--scoreGap " + scoreGap) else ''} \
         ~{if defined(scoreGapNoncan) then ("--scoreGapNoncan " + scoreGapNoncan) else ''} \
         ~{if defined(scoreGapGCAG) then ("--scoreGapGCAG " + scoreGapGCAG) else ''} \
         ~{if defined(scoreGapATAC) then ("--scoreGapATAC " + scoreGapATAC) else ''} \
         ~{if defined(scoreGenomicLengthLog2scale) then ("--scoreGenomicLengthLog2scale " + scoreGenomicLengthLog2scale) else ''} \
         ~{if defined(scoreDelOpen) then ("--scoreDelOpen " + scoreDelOpen) else ''} \
         ~{if defined(scoreDelBase) then ("--scoreDelBase " + scoreDelBase) else ''} \
         ~{if defined(scoreInsOpen) then ("--scoreInsOpen " + scoreInsOpen) else ''} \
         ~{if defined(scoreInsBase) then ("--scoreInsBase " + scoreInsBase) else ''} \
         ~{if defined(scoreStitchSJshift) then ("--scoreStitchSJshift " + scoreStitchSJshift) else ''} \
         ~{if defined(seedSearchStartLmax) then ("--seedSearchStartLmax " + seedSearchStartLmax) else ''} \
         ~{if defined(seedSearchStartLmaxOverLread) then ("--seedSearchStartLmaxOverLread " + seedSearchStartLmaxOverLread) else ''} \
         ~{if defined(seedSearchLmax) then ("--seedSearchLmax " + seedSearchLmax) else ''} \
         ~{if defined(seedMultimapNmax) then ("--seedMultimapNmax " + seedMultimapNmax) else ''} \
         ~{if defined(seedPerReadNmax) then ("--seedPerReadNmax " + seedPerReadNmax) else ''} \
         ~{if defined(seedPerWindowNmax) then ("--seedPerWindowNmax " + seedPerWindowNmax) else ''} \
         ~{if defined(seedNoneLociPerWindow) then ("--seedNoneLociPerWindow " + seedNoneLociPerWindow) else ''} \
         ~{if defined(seedSplitMin) then ("--seedSplitMin " + seedSplitMin) else ''} \
         ~{if defined(alignIntronMin) then ("--alignIntronMin " + alignIntronMin) else ''} \
         ~{if defined(alignIntronMax) then ("--alignIntronMax " + alignIntronMax) else ''} \
         ~{if defined(alignMatesGapMax) then ("--alignMatesGapMax " + alignMatesGapMax) else ''} \
         ~{if defined(alignSJoverhangMin) then ("--alignSJoverhangMin " + alignSJoverhangMin) else ''} \
         ~{if (defined(alignSJstitchMismatchNmax) && length(select_first([alignSJstitchMismatchNmax])) > 0) then "--alignSJstitchMismatchNmax " + sep(" ", select_first([alignSJstitchMismatchNmax])) else ""} \
         ~{if defined(alignSJDBoverhangMin) then ("--alignSJDBoverhangMin " + alignSJDBoverhangMin) else ''} \
         ~{if defined(alignSplicedMateMapLmin) then ("--alignSplicedMateMapLmin " + alignSplicedMateMapLmin) else ''} \
         ~{if defined(alignSplicedMateMapLminOverLmate) then ("--alignSplicedMateMapLminOverLmate " + alignSplicedMateMapLminOverLmate) else ''} \
         ~{if defined(alignWindowsPerReadNmax) then ("--alignWindowsPerReadNmax " + alignWindowsPerReadNmax) else ''} \
         ~{if defined(alignTranscriptsPerWindowNmax) then ("--alignTranscriptsPerWindowNmax " + alignTranscriptsPerWindowNmax) else ''} \
         ~{if defined(alignTranscriptsPerReadNmax) then ("--alignTranscriptsPerReadNmax " + alignTranscriptsPerReadNmax) else ''} \
         ~{if defined(alignEndsType) then ("--alignEndsType '" + alignEndsType + "'") else ""} \
         ~{if defined(alignEndsProtrude) then ("--alignEndsProtrude " + alignEndsProtrude) else ''} \
         ~{if defined(alignSoftClipAtReferenceEnds) then ("--alignSoftClipAtReferenceEnds '" + alignSoftClipAtReferenceEnds + "'") else ""} \
         ~{if defined(alignInsertionFlush) then ("--alignInsertionFlush '" + alignInsertionFlush + "'") else ""} \
         ~{if defined(peOverlapNbasesMin) then ("--peOverlapNbasesMin " + peOverlapNbasesMin) else ''} \
         ~{if defined(peOverlapMMp) then ("--peOverlapMMp " + peOverlapMMp) else ''} \
         ~{if defined(winAnchorMultimapNmax) then ("--winAnchorMultimapNmax " + winAnchorMultimapNmax) else ''} \
         ~{if defined(winBinNbits) then ("--winBinNbits " + winBinNbits) else ''} \
         ~{if defined(winAnchorDistNbins) then ("--winAnchorDistNbins " + winAnchorDistNbins) else ''} \
         ~{if defined(winFlankNbins) then ("--winFlankNbins " + winFlankNbins) else ''} \
         ~{if defined(winReadCoverageRelativeMin) then ("--winReadCoverageRelativeMin " + winReadCoverageRelativeMin) else ''} \
         ~{if defined(winReadCoverageBasesMin) then ("--winReadCoverageBasesMin " + winReadCoverageBasesMin) else ''} \
         ~{if (defined(chimOutType) && length(select_first([chimOutType])) > 0) then "--chimOutType '" + sep("' '", select_first([chimOutType])) + "'" else ""} \
         ~{if defined(chimSegmentMin) then ("--chimSegmentMin " + chimSegmentMin) else ''} \
         ~{if defined(chimScoreMin) then ("--chimScoreMin " + chimScoreMin) else ''} \
         ~{if defined(chimScoreDropMax) then ("--chimScoreDropMax " + chimScoreDropMax) else ''} \
         ~{if defined(chimScoreSeparation) then ("--chimScoreSeparation " + chimScoreSeparation) else ''} \
         ~{if defined(chimScoreJunctionNonGTAG) then ("--chimScoreJunctionNonGTAG " + chimScoreJunctionNonGTAG) else ''} \
         ~{if defined(chimJunctionOverhangMin) then ("--chimJunctionOverhangMin " + chimJunctionOverhangMin) else ''} \
         ~{if defined(chimSegmentReadGapMax) then ("--chimSegmentReadGapMax " + chimSegmentReadGapMax) else ''} \
         ~{if defined(chimFilter) then ("--chimFilter '" + chimFilter + "'") else ""} \
         ~{if defined(chimMainSegmentMultNmax) then ("--chimMainSegmentMultNmax " + chimMainSegmentMultNmax) else ''} \
         ~{if defined(chimMultimapNmax) then ("--chimMultimapNmax " + chimMultimapNmax) else ''} \
         ~{if defined(chimMultimapScoreRange) then ("--chimMultimapScoreRange " + chimMultimapScoreRange) else ''} \
         ~{if defined(chimNonchimScoreDropMin) then ("--chimNonchimScoreDropMin " + chimNonchimScoreDropMin) else ''} \
         ~{if defined(chimOutJunctionFormat) then ("--chimOutJunctionFormat " + chimOutJunctionFormat) else ''} \
         ~{if defined(quantMode) then ("--quantMode '" + quantMode + "'") else ""} \
         ~{if defined(quantTranscriptomeBAMcompression) then ("--quantTranscriptomeBAMcompression " + quantTranscriptomeBAMcompression) else ''} \
         ~{if defined(quantTranscriptomeBan) then ("--quantTranscriptomeBan '" + quantTranscriptomeBan + "'") else ""} \
         ~{if defined(twopassMode) then ("--twopassMode '" + twopassMode + "'") else ""} \
         ~{if defined(twopass1readsN) then ("--twopass1readsN " + twopass1readsN) else ''} \
         ~{if defined(waspOutputMode) then ("--waspOutputMode '" + waspOutputMode + "'") else ""} \
         ~{if defined(soloType) then ("--soloType '" + soloType + "'") else ""} \
         ~{if defined(soloCBwhitelist) then ("--soloCBwhitelist '" + soloCBwhitelist + "'") else ""} \
         ~{if defined(soloCBstart) then ("--soloCBstart " + soloCBstart) else ''} \
         ~{if defined(soloCBlen) then ("--soloCBlen " + soloCBlen) else ''} \
         ~{if defined(soloUMIstart) then ("--soloUMIstart " + soloUMIstart) else ''} \
         ~{if defined(soloUMIlen) then ("--soloUMIlen " + soloUMIlen) else ''} \
         ~{if defined(soloBarcodeReadLength) then ("--soloBarcodeReadLength " + soloBarcodeReadLength) else ''} \
         ~{if (defined(soloCBposition) && length(select_first([soloCBposition])) > 0) then "--soloCBposition '" + sep("' '", select_first([soloCBposition])) + "'" else ""} \
         ~{if defined(soloUMIposition) then ("--soloUMIposition '" + soloUMIposition + "'") else ""} \
         ~{if defined(soloAdapterSequence) then ("--soloAdapterSequence '" + soloAdapterSequence + "'") else ""} \
         ~{if defined(soloAdapterMismatchesNmax) then ("--soloAdapterMismatchesNmax " + soloAdapterMismatchesNmax) else ''} \
         ~{if defined(soloCBmatchWLtype) then ("--soloCBmatchWLtype '" + soloCBmatchWLtype + "'") else ""} \
         ~{if defined(soloStrand) then ("--soloStrand '" + soloStrand + "'") else ""} \
         ~{if defined(soloFeatures) then ("--soloFeatures '" + soloFeatures + "'") else ""} \
         ~{if defined(soloUMIdedup) then ("--soloUMIdedup '" + soloUMIdedup + "'") else ""} \
         ~{if defined(soloUMIfiltering) then ("--soloUMIfiltering '" + soloUMIfiltering + "'") else ""} \
         ~{if defined(soloOutFileNames) then ("--soloOutFileNames '" + soloOutFileNames + "'") else ""} \
         ~{if defined(soloCellFilter) then ("--soloCellFilter '" + soloCellFilter + "'") else ""} \
         --runMode 'genomeGenerate'
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 4, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/star:2.5.3a--0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 32, 4])}G"
       preemptible: 2
     }
     output {
       Directory out = glob(".")[0]
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: STAR Aligner
   doc: |
     Spliced Transcripts Alignment to a Reference © Alexander Dobin, 2009-2019 

     For more details see:

     - https://www.ncbi.nlm.nih.gov/pubmed/23104886
     - https://github.com/alexdobin/STAR
     - https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/star:2.5.3a--0

   inputs:
   - id: parametersFiles
     label: parametersFiles
     doc: '(default: -) none. Can only be defined on the command line.'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --parametersFiles
   - id: sysShell
     label: sysShell
     doc: |-
       (default: -) path to the shell binary, preferably bash, e.g. /bin/bash.
       - ... the default shell is executed, typically /bin/sh. This was reported to fail on some Ubuntu systems - then you need to specify path to bash.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sysShell
   - id: runThreadN
     label: runThreadN
     doc: '(default: 1) number of threads to run STAR'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --runThreadN
       valueFrom: |-
         $([inputs.runtime_cpu, 4, 1].filter(function (inner) { return inner != null })[0])
   - id: runDirPerm
     label: runDirPerm
     doc: |-
       (default: User_RWX) permissions for the directories created at the run-time. 
       - User_RWX ... user-read/write/execute 
       - All_RWX  ... all-read/write/execute (same as chmod 777)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --runDirPerm
   - id: runRNGseed
     label: runRNGseed
     doc: '(default: 777) random number generator seed.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --runRNGseed
   - id: genomeDir
     label: genomeDir
     doc: '(default: GenomeDir/) path to the directory where genome files are stored'
     type:
     - Directory
     - 'null'
     inputBinding:
       prefix: --genomeDir
   - id: outputGenomeDir
     label: outputGenomeDir
     doc: generated for --runMode generateGenome
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genomeDir
   - id: genomeLoad
     label: genomeLoad
     doc: |-
       (default: NoSharedMemory) mode of shared memory usage for the genome files. Only used with --runMode alignReads.
       - LoadAndKeep     ... load genome into shared and keep it in memory after run,
       - LoadAndRemove   ... load genome into shared but remove it after run,
       - LoadAndExit     ... load genome into shared memory and exit, keeping the genome in memory for future runs,
       - Remove:   ... do not map anything, just remove loaded genome from memory,
       - NoSharedMemory  ... do not use shared memory, each job will have its own private copy of the genome
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --genomeLoad
   - id: genomeFastaFiles
     label: genomeFastaFiles
     doc: |-
       (default: -) path(s) to the fasta files with the genome sequences, separated by spaces. These files should be plain text FASTA files, they *cannot* be zipped. Required for the genome generation (--runMode genomeGenerate). Can also be used in the mapping (--runMode alignReads) to add extra (new) sequences to the genome (e.g. spike-ins).
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --genomeFastaFiles
   - id: genomeChainFiles
     label: genomeChainFiles
     doc: |-
       (default: -) chain files for genomic liftover. Only used with --runMode liftOver .
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --genomeChainFiles
   - id: genomeFileSizes
     label: genomeFileSizes
     doc: |-
       (default: 0) genome files exact sizes in bytes. Typically, this should not be defined by the user.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --genomeFileSizes
   - id: genomeConsensusFile
     label: genomeConsensusFile
     doc: |-
       (default: -) VCF file with consensus SNPs (i.e. alternative allele is the major (AF>0.5) allele)
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --genomeConsensusFile
   - id: genomeChrBinNbits
     label: genomeChrBinNbits
     doc: |-
       (default: 18) each chromosome will occupy an integer number of bins. For a genome with large number of contigs, it is recommended to scale this parameter as ``min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)])``.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --genomeChrBinNbits
   - id: genomeSAindexNbases
     label: genomeSAindexNbases
     doc: |-
       (default: 14) length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, the parameter --genomeSAindexNbases must be scaled down to min(14, log2(GenomeLength)/2 - 1).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --genomeSAindexNbases
   - id: genomeSAsparseD
     label: genomeSAsparseD
     doc: |-
       (default: 1) use bigger numbers to decrease needed RAM at the cost of mapping speed reduction
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --genomeSAsparseD
   - id: genomeSuffixLengthMax
     label: genomeSuffixLengthMax
     doc: |-
       (default: -1) maximum length of the suffixes, has to be longer than read length. -1 = infinite.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --genomeSuffixLengthMax
   - id: sjdbFileChrStartEnd
     label: sjdbFileChrStartEnd
     doc: |-
       (default: -) path to the files with genomic coordinates (chr <tab> start <tab> end <tab> strand) for the splice junction introns. Multiple files can be supplied wand will be concatenated.
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --sjdbFileChrStartEnd
   - id: sjdbGTFfile
     label: sjdbGTFfile
     doc: '(default: -) path to the GTF file with annotations'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --sjdbGTFfile
   - id: sjdbGTFchrPrefix
     label: sjdbGTFchrPrefix
     doc: |-
       (default: -) prefix for chromosome names in a GTF file (e.g. 'chr' for using ENSMEBL annotations with UCSC genomes)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFchrPrefix
   - id: sjdbGTFfeatureExon
     label: sjdbGTFfeatureExon
     doc: |-
       (default: exon) feature type in GTF file to be used as exons for building transcripts
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFfeatureExon
   - id: sjdbGTFtagExonParentTranscript
     label: sjdbGTFtagExonParentTranscript
     doc: |-
       (default: transcript_id) GTF attribute name for parent transcript ID (default "transcript_id" works for GTF files)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFtagExonParentTranscript
   - id: sjdbGTFtagExonParentGene
     label: sjdbGTFtagExonParentGene
     doc: |-
       (default: gene_id) GTF attribute name for parent gene ID (default "gene_id" works for GTF files)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFtagExonParentGene
   - id: sjdbGTFtagExonParentGeneName
     label: sjdbGTFtagExonParentGeneName
     doc: '(default: gene_name) GTF attrbute name for parent gene name'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFtagExonParentGeneName
   - id: sjdbGTFtagExonParentGeneType
     label: sjdbGTFtagExonParentGeneType
     doc: '(default: gene_type gene_biotype) GTF attrbute name for parent gene type'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbGTFtagExonParentGeneType
   - id: sjdbOverhang
     label: sjdbOverhang
     doc: |-
       (default: 100) length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --sjdbOverhang
   - id: sjdbScore
     label: sjdbScore
     doc: '(default: 2) extra alignment score for alignmets that cross database junctions'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --sjdbScore
   - id: sjdbInsertSave
     label: sjdbInsertSave
     doc: |-
       (default: Basic) which files to save when sjdb junctions are inserted on the fly at the mapping step
       - Basic ... only small junction / transcript files 
       - All   ... all files including big Genome, SA and SAindex - this will create a complete genome directory
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sjdbInsertSave
   - id: varVCFfile
     label: varVCFfile
     doc: '(default: -) path to the VCF file that contains variation data.'
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --varVCFfile
   - id: inputBAMfile
     label: inputBAMfile
     doc: |-
       (default: -) path to BAM input file, to be used with --runMode inputAlignmentsFromBAM
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --inputBAMfile
   - id: readFilesType
     label: readFilesType
     doc: |
       (default: Fastx) format of input read files
       - Fastx       ... FASTA or FASTQ
       - SAM SE      ... SAM or BAM single-end reads; for BAM use --readFilesCommand samtools view
       - SAM PE      ... SAM or BAM paired-end reads; for BAM use --readFilesCommand samtools view
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readFilesType
   - id: readFilesIn
     label: readFilesIn
     doc: |-
       (default: Read1 Read2) paths to files that contain input read1 (and, if needed,  read2)
     type:
     - type: array
       items: File
     - 'null'
     inputBinding:
       prefix: --readFilesIn
       itemSeparator: ' '
   - id: readFilesPrefix
     label: readFilesPrefix
     doc: |-
       (default: -)   for the read files names, i.e. it will be added in front of the strings in --readFilesIn no prefix
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readFilesPrefix
   - id: readFilesCommand
     label: readFilesCommand
     doc: |-
       (default: -) command line to execute for each of the input file. This command should generate FASTA or FASTQ text and send it to stdout zcat - to uncompress .gz files, bzcat - to uncompress .bz2 files, etc.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readFilesCommand
   - id: readMapNumber
     label: readMapNumber
     doc: '(default: 1) number of reads to map from the beginning of the file map all
       reads'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --readMapNumber
   - id: readMatesLengthsIn
     label: readMatesLengthsIn
     doc: |-
       (default: NotEqual) Equal/NotEqual - lengths of names,sequences,qualities for both mates are the same  / not the same. NotEqual is safe in all situations.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readMatesLengthsIn
   - id: readNameSeparator
     label: readNameSeparator
     doc: |-
       (default: /) character(s) separating the part of the read names that will be trimmed in output (read name after space is always trimmed)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --readNameSeparator
   - id: readQualityScoreBase
     label: readQualityScoreBase
     doc: |-
       (default: 33) number to be subtracted from the ASCII code to get Phred quality score
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --readQualityScoreBase
   - id: clip3pNbases
     label: clip3pNbases
     doc: |-
       (default: 0) number(s) of bases to clip from 3p of each mate. If one value is given, it will be assumed the same for both mates.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --clip3pNbases
   - id: clip5pNbases
     label: clip5pNbases
     doc: |-
       (default: 0) number(s) of bases to clip from 5p of each mate. If one value is given, it will be assumed the same for both mates.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --clip5pNbases
   - id: clip3pAdapterSeq
     label: clip3pAdapterSeq
     doc: |-
       (default: -) adapter sequences to clip from 3p of each mate.  If one value is given, it will be assumed the same for both mates.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --clip3pAdapterSeq
   - id: clip3pAdapterMMp
     label: clip3pAdapterMMp
     doc: |-
       (default: 0.1) max proportion of mismatches for 3p adpater clipping for each mate.  If one value is given, it will be assumed the same for both mates.
     type:
     - double
     - 'null'
     inputBinding:
       prefix: --clip3pAdapterMMp
   - id: clip3pAfterAdapterNbases
     label: clip3pAfterAdapterNbases
     doc: |-
       (default: 0) number of bases to clip from 3p of each mate after the adapter clipping. If one value is given, it will be assumed the same for both mates.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --clip3pAfterAdapterNbases
   - id: limitGenomeGenerateRAM
     label: limitGenomeGenerateRAM
     doc: '(default: 31000000000) maximum available RAM (bytes) for genome generation'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitGenomeGenerateRAM
   - id: limitIObufferSize
     label: limitIObufferSize
     doc: |-
       (default: 150000000) max available buffers size (bytes) for input/output, per thread
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitIObufferSize
   - id: limitOutSAMoneReadBytes
     label: limitOutSAMoneReadBytes
     doc: '(default: 100000) >(2*(LengthMate1+LengthMate2+100)*outFilterMultimapNmax'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitOutSAMoneReadBytes
   - id: limitOutSJoneRead
     label: limitOutSJoneRead
     doc: |-
       (default: 1000) max number of junctions for one read (including all multi-mappers)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitOutSJoneRead
   - id: limitOutSJcollapsed
     label: limitOutSJcollapsed
     doc: '(default: 1000000) max number of collapsed junctions'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitOutSJcollapsed
   - id: limitBAMsortRAM
     label: limitBAMsortRAM
     doc: |-
       (default: 0) maximum available RAM (bytes) for sorting BAM. If =0, it will be set to the genome index size. 0 value can only be used with --genomeLoad NoSharedMemory option.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitBAMsortRAM
   - id: limitSjdbInsertNsj
     label: limitSjdbInsertNsj
     doc: |-
       (default: 1000000) maximum number of junction to be inserted to the genome on the fly at the mapping stage, including those from annotations and those detected in the 1st step of the 2-pass run
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitSjdbInsertNsj
   - id: limitNreadsSoft
     label: limitNreadsSoft
     doc: '(default: 1) soft limit on the number of reads'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --limitNreadsSoft
   - id: outFileNamePrefix
     label: outFileNamePrefix
     doc: |-
       (default: ./) output files name prefix (including full or relative path). Can only be defined on the command line.
     type: string
     default: ./
     inputBinding:
       prefix: --outFileNamePrefix
   - id: outTmpDir
     label: outTmpDir
     doc: |-
       (default: -) path to a directory that will be used as temporary by STAR. All contents of this directory will be removed!     - the temp directory will default to outFileNamePrefix_STARtmp
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outTmpDir
   - id: outTmpKeep
     label: outTmpKeep
     doc: |-
       (default: None) whether to keep the tempporary files after STAR runs is finished None ... remove all temporary files All .. keep all files
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outTmpKeep
   - id: outStd
     label: outStd
     doc: |-
       (default: Log) which output will be directed to stdout (standard out) Log     ... log messages SAM        ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out BAM_Unsorted     ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate BAM_Quant        ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM       
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outStd
   - id: outReadsUnmapped
     label: outReadsUnmapped
     doc: |-
       (default: None) which output will be directed to stdout (standard out) [Log ... log messages SAM ... alignments in SAM format (which normally are output to Aligned.out.sam file), normal standard output will go into Log.std.out BAM_Unsorted           ... alignments in BAM format, unsorted. Requires --outSAMtype BAM Unsorted BAM_SortedByCoordinate ... alignments in BAM format, unsorted. Requires --outSAMtype BAM SortedByCoordinate BAM_Quant ... alignments to transcriptome in BAM format, unsorted. Requires --quantMode TranscriptomeSAM]
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outReadsUnmapped
   - id: outQSconversionAdd
     label: outQSconversionAdd
     doc: |-
       (default: 0) add this number to the quality score (e.g. to convert from Illumina to Sanger, use -31)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outQSconversionAdd
   - id: outMultimapperOrder
     label: outMultimapperOrder
     doc: '(default: Old_2.4) order of multimapping alignments in the output files Old_2.4'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outMultimapperOrder
   - id: outSAMtype
     label: outSAMtype
     doc: |-
       (default: SAM) ... quasi-random order used before 2.5.0 Random ... random order of alignments for each multi-mapper. Read mates (pairs) are always adjacent, all alignment for each read stay together. This option will become default in the future releases. ... standard unsorted SortedByCoordinate ... sorted by coordinate. This option will allocate extra memory for sorting which can be specified by --limitBAMsortRAM.
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --outSAMtype
       itemSeparator: ' '
   - id: outSAMmode
     label: outSAMmode
     doc: |-
       (default: Full) mode of SAM output None ... no SAM output Full ... full SAM output NoQS ... full SAM but without quality scores ... no attributes Standard    ... NH HI AS nM All   ... NH HI AS nM NM MD jM jI MC ch vA    ... variant allele vG    ... genomic coordiante of the variant overlapped by the read vW    ... 0/1 - alignment does not pass / passes WASP filtering. Requires --waspOutputMode SAMtag STARsolo: CR CY UR UY ... sequences and quality scores of cell barcodes and UMIs for the solo* demultiplexing CB UB       ... error-corrected cell barcodes and UMIs for solo* demultiplexing. Requires --outSAMtype BAM SortedByCoordinate. sM    ... assessment of CB and UMI sS    ... sequence of the entire barcode (CB,UMI,adapter...) sQ    ... quality of the entire barcode Unsupported/undocumented: rB    ... alignment block read/genomic coordinates vR    ... read coordinate of the variant
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMmode
   - id: outSAMstrandField
     label: outSAMstrandField
     doc: '(default: None) Cufflinks-like strand field flag None'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMstrandField
   - id: outSAMattributes
     label: outSAMattributes
     doc: |-
       (default: Standard) a string of desired SAM attributes, in the order desired for the output SAM NH HI AS nM NM MD jM jI XS MC ch ... any combination in any order None
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMattributes
   - id: outSAMattrIHstart
     label: outSAMattrIHstart
     doc: |-
       (default: 1) start value for the IH attribute. 0 may be required by some downstream software, such as Cufflinks or StringTie.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMattrIHstart
   - id: outSAMunmapped
     label: outSAMunmapped
     doc: |-
       (default: None) output of unmapped reads in the SAM format 1st word: None   ... no output Within ... output unmapped reads within the main SAM file (i.e. Aligned.out.sam) 2nd word: KeepPairs ... record unmapped mate for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate. Only affects multi-mapping reads.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMunmapped
   - id: outSAMorder
     label: outSAMorder
     doc: |-
       (default: Paired) type of sorting for the SAM output one mate after the other for all paired alignments one mate after the other for all paired alignments, the order is kept the same as in the input FASTQ files
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMorder
   - id: outSAMprimaryFlag
     label: outSAMprimaryFlag
     doc: |-
       (default: OneBestScore) which alignments are considered primary - all others will be marked with 0x100 bit in the FLAG OneBestScore ... only one alignment with the best score is primary AllBestScore ... all alignments with the best score are primary
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMprimaryFlag
   - id: outSAMreadID
     label: outSAMreadID
     doc: |-
       (default: Standard) read ID record type Standard ... first word (until space) from the FASTx read ID line, removing /1,/2 from the end Number   ... read number (index) in the FASTx file
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMreadID
   - id: outSAMmapqUnique
     label: outSAMmapqUnique
     doc: '(default: 255) the MAPQ value for unique mappers'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMmapqUnique
   - id: outSAMflagOR
     label: outSAMflagOR
     doc: |-
       (default: 0) sam FLAG will be bitwise OR'd with this value, i.e. FLAG=FLAG | outSAMflagOR. This is applied after all flags have been set by STAR, and after outSAMflagAND. Can be used to set specific bits that are not set otherwise.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMflagOR
   - id: outSAMflagAND
     label: outSAMflagAND
     doc: |-
       (default: 65535) sam FLAG will be bitwise AND'd with this value, i.e. FLAG=FLAG & outSAMflagOR. This is applied after all flags have been set by STAR, but before outSAMflagOR. Can be used to unset specific bits that are not set otherwise.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMflagAND
   - id: outSAMattrRGline
     label: outSAMattrRGline
     doc: |-
       (default: -) SAM/BAM read group line. The first word contains the read group identifier and must start with "ID:", e.g. --outSAMattrRGline ID:xxx CN:yy "DS:z z z".     xxx will be added as RG tag to each output alignment. Any spaces in the tag values have to be double quoted.     Comma separated RG lines correspons to different (comma separated) input files in --readFilesIn. Commas have to be surrounded by spaces, e.g.     --outSAMattrRGline ID:xxx , ID:zzz "DS:z z" , ID:yyy DS:yyyy
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMattrRGline
   - id: outSAMheaderHD
     label: outSAMheaderHD
     doc: '(default: -) @HD (header) line of the SAM header'
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --outSAMheaderHD
   - id: outSAMheaderPG
     label: outSAMheaderPG
     doc: '(default: -) extra @PG (software) line of the SAM header (in addition to STAR)'
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --outSAMheaderPG
   - id: outSAMheaderCommentFile
     label: outSAMheaderCommentFile
     doc: '(default: -) path to the file with @CO (comment) lines of the SAM header'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMheaderCommentFile
   - id: outSAMfilter
     label: outSAMfilter
     doc: |-
       (default: None) filter the output into main SAM/BAM files KeepOnlyAddedReferences ... only keep the reads for which all alignments are to the extra reference sequences added with --genomeFastaFiles at the mapping stage. KeepAllAddedReferences ...  keep all alignments to the extra reference sequences added with --genomeFastaFiles at the mapping stage.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSAMfilter
   - id: outSAMmultNmax
     label: outSAMmultNmax
     doc: |-
       (default: 1) max number of multiple alignments for a read that will be output to the SAM/BAM files. -1 ... all alignments (up to --outFilterMultimapNmax) will be output
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMmultNmax
   - id: outSAMtlen
     label: outSAMtlen
     doc: |-
       (default: 1) calculation method for the TLEN field in the SAM/BAM files 1 ... leftmost base of the (+)strand mate to rightmost base of the (-)mate. (+)sign for the (+)strand mate 2 ... leftmost base of any mate to rightmost base of any mate. (+)sign for the mate with the leftmost base. This is different from 1 for overlapping mates with protruding ends
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSAMtlen
   - id: outBAMcompression
     label: outBAMcompression
     doc: |-
       (default: 1) -1 to 10  BAM compression level, -1=default compression (6?), 0=no compression, 10=maximum compression
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outBAMcompression
   - id: outBAMsortingThreadN
     label: outBAMsortingThreadN
     doc: |-
       (default: 0) number of threads for BAM sorting. 0 will default to min(6,--runThreadN).
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outBAMsortingThreadN
   - id: outBAMsortingBinsN
     label: outBAMsortingBinsN
     doc: '(default: 50) number of genome bins fo coordinate-sorting'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outBAMsortingBinsN
   - id: bamRemoveDuplicatesType
     label: bamRemoveDuplicatesType
     doc: |-
       (default: -) mark duplicates in the BAM file, for now only works with (i) sorted BAM fed with inputBAMfile, and (ii) for paired-end alignments only -
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --bamRemoveDuplicatesType
   - id: bamRemoveDuplicatesMate2basesN
     label: bamRemoveDuplicatesMate2basesN
     doc: |-
       (default: 0) number of bases from the 5' of mate 2 to use in collapsing (e.g. for RAMPAGE)
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --bamRemoveDuplicatesMate2basesN
   - id: outWigType
     label: outWigType
     doc: |-
       (default: None) --outSAMtype BAM SortedByCoordinate .     1st word:     None       ... no signal output     bedGraph   ... bedGraph format     wiggle     ... wiggle format     2nd word:     read1_5p   ... signal from only 5' of the 1st read, useful for CAGE/RAMPAGE etc     read2      ... signal from only 2nd read
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outWigType
   - id: outWigStrand
     label: outWigStrand
     doc: |-
       (default: Stranded) strandedness of wiggle/bedGraph output     Stranded   ...  separate strands, str1 and str2     Unstranded ...  collapsed strands
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outWigStrand
   - id: outWigReferencesPrefix
     label: outWigReferencesPrefix
     doc: |-
       (default: -) prefix matching reference names to include in the output wiggle file, e.g. "chr", default "-" - include all references
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outWigReferencesPrefix
   - id: outWigNorm
     label: outWigNorm
     doc: |-
       (default: RPM) type of normalization for the signal RPM    ... reads per million of mapped reads None   ... no normalization, "raw" counts
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outWigNorm
   - id: outFilterType
     label: outFilterType
     doc: |-
       (default: Normal) type of filtering Normal  ... standard filtering using only current alignment BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outFilterType
   - id: outFilterMultimapScoreRange
     label: outFilterMultimapScoreRange
     doc: '(default: 1) the score range below the maximum score for multimapping alignments'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outFilterMultimapScoreRange
   - id: outFilterMultimapNmax
     label: outFilterMultimapNmax
     doc: |-
       (default: 10) maximum number of loci the read is allowed to map to. Alignments (all of them) will be output only if the read maps to no more loci than this value.  Otherwise no alignments will be output, and the read will be counted as "mapped to too many loci" in the Log.final.out .
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outFilterMultimapNmax
   - id: outFilterMismatchNmax
     label: outFilterMismatchNmax
     doc: |-
       (default: 10) alignment will be output only if it has no more mismatches than this value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outFilterMismatchNmax
   - id: outFilterMismatchNoverLmax
     label: outFilterMismatchNoverLmax
     doc: |-
       (default: 0.3) alignment will be output only if its ratio of mismatches to *mapped* length is less than or equal to this value.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --outFilterMismatchNoverLmax
   - id: outFilterMismatchNoverReadLmax
     label: outFilterMismatchNoverReadLmax
     doc: |-
       (default: 1) alignment will be output only if its ratio of mismatches to *read* length is less than or equal to this value.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --outFilterMismatchNoverReadLmax
   - id: outFilterScoreMin
     label: outFilterScoreMin
     doc: |-
       (default: 0) alignment will be output only if its score is higher than or equal to this value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outFilterScoreMin
   - id: outFilterScoreMinOverLread
     label: outFilterScoreMinOverLread
     doc: |-
       (default: 0.66) same as outFilterScoreMin, but  normalized to read length (sum of mates' lengths for paired-end reads)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --outFilterScoreMinOverLread
   - id: outFilterMatchNmin
     label: outFilterMatchNmin
     doc: |-
       (default: 0) alignment will be output only if the number of matched bases is higher than or equal to this value.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outFilterMatchNmin
   - id: outFilterMatchNminOverLread
     label: outFilterMatchNminOverLread
     doc: |-
       (default: 0.66) sam as outFilterMatchNmin, but normalized to the read length (sum of mates' lengths for paired-end reads).
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --outFilterMatchNminOverLread
   - id: outFilterIntronMotifs
     label: outFilterIntronMotifs
     doc: '(default: None) filter alignment using their motifs None'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outFilterIntronMotifs
   - id: outFilterIntronStrands
     label: outFilterIntronStrands
     doc: |-
       (default: RemoveInconsistentStrands) filter alignments RemoveInconsistentStrands      ... remove alignments that have junctions with inconsistent strands None
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outFilterIntronStrands
   - id: outSJfilterReads
     label: outSJfilterReads
     doc: |-
       (default: All) which reads to consider for collapsed splice junctions output all reads, unique- and multi-mappers uniquely mapping reads only
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --outSJfilterReads
   - id: outSJfilterOverhangMin
     label: outSJfilterOverhangMin
     doc: |-
       (default: 30 12 12 12) minimum overhang length for splice junctions on both sides for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif does not apply to annotated junctions
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSJfilterOverhangMin
   - id: outSJfilterCountUniqueMin
     label: outSJfilterCountUniqueMin
     doc: |-
       (default: 3 1 1 1) minimum uniquely mapping read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied does not apply to annotated junctions
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSJfilterCountUniqueMin
   - id: outSJfilterCountTotalMin
     label: outSJfilterCountTotalMin
     doc: |-
       (default: 3 1 1 1) minimum total (multi-mapping+unique) read count per junction for: (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif. -1 means no output for that motif Junctions are output if one of outSJfilterCountUniqueMin OR outSJfilterCountTotalMin conditions are satisfied does not apply to annotated junctions
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSJfilterCountTotalMin
   - id: outSJfilterDistToOtherSJmin
     label: outSJfilterDistToOtherSJmin
     doc: |-
       (default: 10 0 5 10) minimum allowed distance to other junctions' donor/acceptor does not apply to annotated junctions
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSJfilterDistToOtherSJmin
   - id: outSJfilterIntronMaxVsReadN
     label: outSJfilterIntronMaxVsReadN
     doc: |-
       (default: 50000 100000 200000) maximum gap allowed for junctions supported by 1,2,3,,,N reads <=200000. by >=4 reads any gap <=alignIntronMax does not apply to annotated junctions
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --outSJfilterIntronMaxVsReadN
   - id: scoreGap
     label: scoreGap
     doc: '(default: 0) splice junction penalty (independent on intron motif)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreGap
   - id: scoreGapNoncan
     label: scoreGapNoncan
     doc: '(default: 8) non-canonical junction penalty (in addition to scoreGap)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreGapNoncan
   - id: scoreGapGCAG
     label: scoreGapGCAG
     doc: '(default: 4) GC/AG and CT/GC junction penalty (in addition to scoreGap)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreGapGCAG
   - id: scoreGapATAC
     label: scoreGapATAC
     doc: '(default: 8) AT/AC  and GT/AT junction penalty  (in addition to scoreGap)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreGapATAC
   - id: scoreGenomicLengthLog2scale
     label: scoreGenomicLengthLog2scale
     doc: '(default: -0.25) scoreGenomicLengthLog2scale*log2(genomicLength)'
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --scoreGenomicLengthLog2scale
   - id: scoreDelOpen
     label: scoreDelOpen
     doc: '(default: 2) deletion open penalty'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreDelOpen
   - id: scoreDelBase
     label: scoreDelBase
     doc: '(default: 2) deletion extension penalty per base (in addition to scoreDelOpen)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreDelBase
   - id: scoreInsOpen
     label: scoreInsOpen
     doc: '(default: 2) insertion open penalty'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreInsOpen
   - id: scoreInsBase
     label: scoreInsBase
     doc: '(default: 2) insertion extension penalty per base (in addition to scoreInsOpen)'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreInsBase
   - id: scoreStitchSJshift
     label: scoreStitchSJshift
     doc: |-
       (default: 1) maximum score reduction while searching for SJ boundaries inthe stitching step
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --scoreStitchSJshift
   - id: seedSearchStartLmax
     label: seedSearchStartLmax
     doc: |-
       (default: 50) defines the search start point through the read - the read is split into pieces no longer than this value
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedSearchStartLmax
   - id: seedSearchStartLmaxOverLread
     label: seedSearchStartLmaxOverLread
     doc: |-
       (default: 1) seedSearchStartLmax normalized to read length (sum of mates' lengths for paired-end reads)
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --seedSearchStartLmaxOverLread
   - id: seedSearchLmax
     label: seedSearchLmax
     doc: |-
       (default: 0) defines the maximum length of the seeds, if =0 max seed lengthis infinite
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedSearchLmax
   - id: seedMultimapNmax
     label: seedMultimapNmax
     doc: |-
       (default: 10000) only pieces that map fewer than this value are utilized in the stitching procedure
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedMultimapNmax
   - id: seedPerReadNmax
     label: seedPerReadNmax
     doc: '(default: 1000) max number of seeds per read'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedPerReadNmax
   - id: seedPerWindowNmax
     label: seedPerWindowNmax
     doc: '(default: 50) max number of seeds per window'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedPerWindowNmax
   - id: seedNoneLociPerWindow
     label: seedNoneLociPerWindow
     doc: '(default: 10) max number of one seed loci per window'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedNoneLociPerWindow
   - id: seedSplitMin
     label: seedSplitMin
     doc: '(default: 12) min length of the seed sequences split by Ns or mate gap'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --seedSplitMin
   - id: alignIntronMin
     label: alignIntronMin
     doc: |-
       (default: 21) genomic gap is considered intron if its length>=alignIntronMin, otherwise it is considered Deletion
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignIntronMin
   - id: alignIntronMax
     label: alignIntronMax
     doc: |-
       (default: 0) maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignIntronMax
   - id: alignMatesGapMax
     label: alignMatesGapMax
     doc: |-
       (default: 0) maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignMatesGapMax
   - id: alignSJoverhangMin
     label: alignSJoverhangMin
     doc: '(default: 5) minimum overhang (i.e. block size) for spliced alignments'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignSJoverhangMin
   - id: alignSJstitchMismatchNmax
     label: alignSJstitchMismatchNmax
     doc: |-
       (default: 0 -1 0 0) maximum number of mismatches for stitching of the splice junctions (-1: no limit).     (1) non-canonical motifs, (2) GT/AG and CT/AC motif, (3) GC/AG and CT/GC motif, (4) AT/AC and GT/AT motif.
     type:
     - type: array
       items: int
     - 'null'
     inputBinding:
       prefix: --alignSJstitchMismatchNmax
   - id: alignSJDBoverhangMin
     label: alignSJDBoverhangMin
     doc: |-
       (default: 3) minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignSJDBoverhangMin
   - id: alignSplicedMateMapLmin
     label: alignSplicedMateMapLmin
     doc: '(default: 0) minimum mapped length for a read mate that is spliced'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignSplicedMateMapLmin
   - id: alignSplicedMateMapLminOverLmate
     label: alignSplicedMateMapLminOverLmate
     doc: '(default: 0.66) alignSplicedMateMapLmin normalized to mate length'
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --alignSplicedMateMapLminOverLmate
   - id: alignWindowsPerReadNmax
     label: alignWindowsPerReadNmax
     doc: '(default: 10000) max number of windows per read'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignWindowsPerReadNmax
   - id: alignTranscriptsPerWindowNmax
     label: alignTranscriptsPerWindowNmax
     doc: '(default: 100) max number of transcripts per window'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignTranscriptsPerWindowNmax
   - id: alignTranscriptsPerReadNmax
     label: alignTranscriptsPerReadNmax
     doc: '(default: 10000) max number of different alignments per read to consider'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignTranscriptsPerReadNmax
   - id: alignEndsType
     label: alignEndsType
     doc: '(default: Local) type of read ends alignment Local'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --alignEndsType
   - id: alignEndsProtrude
     label: alignEndsProtrude
     doc: |-
       (default: 0 ConcordantPair) allow protrusion of alignment ends, i.e. start (end) of the +strand mate downstream of the start (end) of the -strand mate maximum number of protrusion bases allowed string:     ConcordantPair ... report alignments with non-zero protrusion as concordant pairs     DiscordantPair ... report alignments with non-zero protrusion as discordant pairs
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --alignEndsProtrude
   - id: alignSoftClipAtReferenceEnds
     label: alignSoftClipAtReferenceEnds
     doc: |-
       (default: Yes) allow the soft-clipping of the alignments past the end of the chromosomes Yes ... allow No  ... prohibit, useful for compatibility with Cufflinks
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --alignSoftClipAtReferenceEnds
   - id: alignInsertionFlush
     label: alignInsertionFlush
     doc: |-
       (default: None) how to flush ambiguous insertion positions None    ... insertions are not flushed Right   ... insertions are flushed to the right
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --alignInsertionFlush
   - id: peOverlapNbasesMin
     label: peOverlapNbasesMin
     doc: |-
       (default: 0) minimum number of overlap bases to trigger mates merging and realignment
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --peOverlapNbasesMin
   - id: peOverlapMMp
     label: peOverlapMMp
     doc: '(default: 0.01) maximum proportion of mismatched bases in the overlap area'
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --peOverlapMMp
   - id: winAnchorMultimapNmax
     label: winAnchorMultimapNmax
     doc: '(default: 50) max number of loci anchors are allowed to map to'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --winAnchorMultimapNmax
   - id: winBinNbits
     label: winBinNbits
     doc: |-
       (default: 16) =LOG2(winBin), where winBin is the size of the bin for the windows/clustering, each window will occupy an integer number of bins.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --winBinNbits
   - id: winAnchorDistNbins
     label: winAnchorDistNbins
     doc: |-
       (default: 9) max number of bins between two anchors that allows aggregation of anchors into one window
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --winAnchorDistNbins
   - id: winFlankNbins
     label: winFlankNbins
     doc: |-
       (default: 4) log2(winFlank), where win Flank is the size of the left and right flanking regions for each window
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --winFlankNbins
   - id: winReadCoverageRelativeMin
     label: winReadCoverageRelativeMin
     doc: |-
       (default: 0.5) minimum relative coverage of the read sequence by the seeds in a window, for STARlong algorithm only.
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --winReadCoverageRelativeMin
   - id: winReadCoverageBasesMin
     label: winReadCoverageBasesMin
     doc: |-
       (default: 0) minimum number of bases covered by the seeds in a window , for STARlong algorithm only.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --winReadCoverageBasesMin
   - id: chimOutType
     label: chimOutType
     doc: |-
       (default: Junctions) type of chimeric output     Junctions       ... Chimeric.out.junction     SeparateSAMold  ... output old SAM into separate Chimeric.out.sam file     WithinBAM       ... output into main aligned BAM files (Aligned.*.bam)     WithinBAM HardClip  ... (default) hard-clipping in the CIGAR for supplemental chimeric alignments (defaultif no 2nd word is present)     WithinBAM SoftClip  ... soft-clipping in the CIGAR for supplemental chimeric alignments
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --chimOutType
   - id: chimSegmentMin
     label: chimSegmentMin
     doc: |-
       (default: 0) minimum length of chimeric segment length, if ==0, no chimeric output
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimSegmentMin
   - id: chimScoreMin
     label: chimScoreMin
     doc: '(default: 0) minimum total (summed) score of the chimeric segments'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimScoreMin
   - id: chimScoreDropMax
     label: chimScoreDropMax
     doc: |-
       (default: 20) max drop (difference) of chimeric score (the sum of scores of all chimeric segments) from the read length
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimScoreDropMax
   - id: chimScoreSeparation
     label: chimScoreSeparation
     doc: |-
       (default: 10) minimum difference (separation) between the best chimeric score and the next one
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimScoreSeparation
   - id: chimScoreJunctionNonGTAG
     label: chimScoreJunctionNonGTAG
     doc: '(default: -1) penalty for a non-GT/AG chimeric junction'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimScoreJunctionNonGTAG
   - id: chimJunctionOverhangMin
     label: chimJunctionOverhangMin
     doc: '(default: 20) minimum overhang for a chimeric junction'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimJunctionOverhangMin
   - id: chimSegmentReadGapMax
     label: chimSegmentReadGapMax
     doc: '(default: 0) maximum gap in the read sequence between chimeric segments'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimSegmentReadGapMax
   - id: chimFilter
     label: chimFilter
     doc: |-
       (default: banGenomicN) different filters for chimeric alignments     None ... no filtering     banGenomicN ... Ns are not allowed in the genome sequence around the chimeric junction
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --chimFilter
   - id: chimMainSegmentMultNmax
     label: chimMainSegmentMultNmax
     doc: |-
       (default: 10) maximum number of multi-alignments for the main chimeric segment. =1 will prohibit multimapping main segments.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimMainSegmentMultNmax
   - id: chimMultimapNmax
     label: chimMultimapNmax
     doc: |-
       (default: 0) maximum number of chimeric multi-alignments 0 ... use the old scheme for chimeric detection which only considered unique alignments
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimMultimapNmax
   - id: chimMultimapScoreRange
     label: chimMultimapScoreRange
     doc: |-
       (default: 1) the score range for multi-mapping chimeras below the best chimeric score. Only works with --chimMultimapNmax > 1
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimMultimapScoreRange
   - id: chimNonchimScoreDropMin
     label: chimNonchimScoreDropMin
     doc: |-
       (default: 20) to trigger chimeric detection, the drop in the best non-chimeric alignment score with respect to the read length has to be greater than this value ... none     TranscriptomeSAM ... output SAM/BAM alignments to transcriptome into a separate file     GeneCounts       ... count reads per gene
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimNonchimScoreDropMin
   - id: chimOutJunctionFormat
     label: chimOutJunctionFormat
     doc: |-
       (default: 0) formatting type for the Chimeric.out.junction file 0 ... no comment lines/headers total, unique, multi
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --chimOutJunctionFormat
   - id: quantMode
     label: quantMode
     doc: |-
       (default: -) types of quantification requested     -        ... prohibit single-end alignments
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --quantMode
   - id: quantTranscriptomeBAMcompression
     label: quantTranscriptomeBAMcompression
     doc: |-
       (default: 1 1) -2 to 10  transcriptome BAM compression level     -2  ... no BAM output     -1  ... default compression (6?)      0  ... no compression      10 ... maximum compression ... 1-pass mapping     Basic       ... basic 2-pass mapping, with all 1st pass junctions inserted into the genome indices on the fly
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --quantTranscriptomeBAMcompression
   - id: quantTranscriptomeBan
     label: quantTranscriptomeBan
     doc: |-
       (default: IndelSoftclipSingleend) prohibit various alignment type     IndelSoftclipSingleend  ... prohibit indels, soft clipping and single-end alignments - compatible with RSEM     Singleend
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --quantTranscriptomeBan
   - id: twopassMode
     label: twopassMode
     doc: '(default: None) 2-pass mapping mode.     None'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --twopassMode
   - id: twopass1readsN
     label: twopass1readsN
     doc: |-
       (default: 1) number of reads to process for the 1st step. Use very large number (or default -1) to map all reads in the first step.
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --twopass1readsN
   - id: waspOutputMode
     label: waspOutputMode
     doc: |-
       (default: None) Nature Methods 12, 1061–1063 (2015), https://www.nature.com/articles/nmeth.3582 .     SAMtag      ... add WASP tags to the alignments that pass WASP filtering
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --waspOutputMode
   - id: soloType
     label: soloType
     doc: |-
       (default: None) type of single-cell RNA-seq     CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium     CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloType
   - id: soloCBwhitelist
     label: soloCBwhitelist
     doc: |-
       (default: -) file(s) with whitelist(s) of cell barcodes. Only one file allowed with 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloCBwhitelist
   - id: soloCBstart
     label: soloCBstart
     doc: '(default: 1) cell barcode start base'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloCBstart
   - id: soloCBlen
     label: soloCBlen
     doc: '(default: 16) cell barcode length'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloCBlen
   - id: soloUMIstart
     label: soloUMIstart
     doc: '(default: 17) UMI start base'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloUMIstart
   - id: soloUMIlen
     label: soloUMIlen
     doc: '(default: 10) UMI length'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloUMIlen
   - id: soloBarcodeReadLength
     label: soloBarcodeReadLength
     doc: |-
       (default: 1) length of the barcode read     1   ... equal to sum of soloCBlen+soloUMIlen     0   ... not defined, do not check
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloBarcodeReadLength
   - id: soloCBposition
     label: soloCBposition
     doc: |-
       (default: -) position of Cell Barcode(s) on the barcode read.     Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2. startAnchor_startDistance_endAnchor_endDistance adapter end     start(end)Distance is the distance from the CB start(end) to the Anchor base     String for different barcodes are separated by space. inDrop (Zilionis et al, Nat. Protocols, 2017):     --soloCBposition  0_0_2_-1  3_1_3_8
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --soloCBposition
   - id: soloUMIposition
     label: soloUMIposition
     doc: |-
       (default: -) position of the UMI on the barcode read, same as soloCBposition inDrop (Zilionis et al, Nat. Protocols, 2017):     --soloCBposition  3_9_3_14
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloUMIposition
   - id: soloAdapterSequence
     label: soloAdapterSequence
     doc: |-
       (default: -) adapter sequence to anchor barcodes.    ... only exact matches allowed     1MM         ... only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.     1MM_multi         ... multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches.  Allowed CBs have to have at least one read with exact match. Similar to CellRanger 2.2.0     1MM_multi_pseudocounts  ... same as 1MM_Multi, but pseudocounts of 1 are added to all whitelist barcodes. Similar to CellRanger 3.x.x 
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloAdapterSequence
   - id: soloAdapterMismatchesNmax
     label: soloAdapterMismatchesNmax
     doc: '(default: 1) maximum number of mismatches allowed in adapter sequence.'
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --soloAdapterMismatchesNmax
   - id: soloCBmatchWLtype
     label: soloCBmatchWLtype
     doc: '(default: 1MM_multi) matching the Cell Barcodes to the WhiteList     Exact'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloCBmatchWLtype
   - id: soloStrand
     label: soloStrand
     doc: |-
       (default: Forward) strandedness of the solo libraries:     Unstranded  ... no strand information     Forward     ... read strand same as the original RNA molecule     Reverse     ... read strand opposite to the original RNA molecule .. all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)     1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).     Exact       ... only exactly matching UMIs are collapsed
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloStrand
   - id: soloFeatures
     label: soloFeatures
     doc: |-
       (default: Gene) genomic features for which the UMI counts per Cell Barcode are collected reads match the gene transcript reported in SJ.out.tab count all reads overlapping genes' exons and introns     Transcript3p   ... quantification of transcript for 3' protocols
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloFeatures
   - id: soloUMIdedup
     label: soloUMIdedup
     doc: '(default: 1MM_All) type of UMI deduplication (collapsing) algorithm     1MM_All'
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloUMIdedup
   - id: soloUMIfiltering
     label: soloUMIfiltering
     doc: |-
       (default: -) type of UMI filtering remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)     MultiGeneUMI    ... remove lower-count UMIs that map to more than one gene (introduced in CellRanger 3.x.x)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloUMIfiltering
   - id: soloOutFileNames
     label: soloOutFileNames
     doc: |-
       (default: Solo.out/  features.tsv barcodes.tsv matrix.mtx) file names for STARsolo output:     file_name_prefix   gene_names   barcode_sequences   cell_feature_count_matrix
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloOutFileNames
   - id: soloCellFilter
     label: soloCellFilter
     doc: |-
       (default: CellRanger2.2 3000 0.99 10) ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)     1MM_Directional     ... follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).     Exact       ... only exactly matching UMIs are collapsed
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --soloCellFilter

   outputs:
   - id: out
     label: out
     type: Directory
     outputBinding:
       glob: .
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: STAR
   arguments:
   - prefix: --runMode
     position: 0
     valueFrom: genomeGenerate

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: star_genomeGenerate


