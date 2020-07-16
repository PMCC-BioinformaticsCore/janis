:orphan:

Vep (Cache)
=================

*0 contributors · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.ensembl.vep.v98_3 import VepCache_98_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vep_step",
           VepCache_98_3(
               inputFile=None,
               vcf=None,
           )
       )
       wf.output("std", source=vep_step.std)
       wf.output("out", source=vep_step.out)
       wf.output("stats", source=vep_step.stats)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vep:

.. code-block:: bash

   # user inputs
   janis inputs vep > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputFile: inputFile.vcf.gz




5. Run vep with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vep





Information
------------


:ID: ``vep``
:URL: *No URL to the documentation was provided*
:Versions: 98.3
:Container: ensemblorg/ensembl-vep:release_98.3
:Authors: 
:Citations: None
:Created: None
:Updated: None



Outputs
-----------

======  ============  ===============
name    type          documentation
======  ============  ===============
std     stdout<File>
out     File
stats   File
======  ============  ===============



Additional configuration (inputs)
---------------------------------

====================  ==========================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                        prefix                     position    documentation
====================  ==========================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
inputFile             CompressedVCF               --input_file                           Input file name. Can use compressed file (gzipped).
vcf                   Boolean                     --vcf                                  Writes output in VCF format. Consequences are added in the INFO field of the VCF file, using the key "CSQ". Data fields are encoded separated by "|"; the order of fields is written in the VCF header. Output fields in the "CSQ" INFO field can be selected by using --fields. If the input format was VCF, the file will remain unchanged save for the addition of the CSQ field (unless using any filtering). Custom data added with --custom are added as separate fields, using the key specified for each data file. Commas in fields are replaced with ampersands (&) to preserve VCF format.
outputFilename        Optional<Filename>          --output_file                          (-o) Output file name. Results can write to STDOUT by specifying  as the output file name - this will force quiet mode. Default = "variant_effect_output.txt"
help                  Optional<Boolean>           --help                                 Display help message and quit
quiet                 Optional<Boolean>           --quiet                                (-q) Suppress warning messages.Not used by default
verbose               Optional<Boolean>           --verbose                              (-v) Print out a bit more information while running. Not used by default
config                Optional<File>              --config                               Load configuration options from a config file. The config file should consist of whitespace-separated pairs of option names and settings e.g.:

                                                                                                     output_file   my_output.txt
                                                                                                     species       mus_musculus
                                                                                                     format        vcf
                                                                                                     host          useastdb.ensembl.org

                                                                                                     A config file can also be implicitly read; save the file as $HOME/.vep/vep.ini (or equivalent directory if
                                                                                                     using --dir). Any options in this file will be overridden by those specified in a config file using --config,
                                                                                                     and in turn by any options specified on the command line. You can create a quick version file of this by
                                                                                                     setting the flags as normal and running VEP in verbose (-v) mode. This will output lines that can be copied
                                                                                                     to a config file that can be loaded in on the next run using --config. Not used by default
everything            Optional<Boolean>           --everything                           (-e) Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --variant_class, --mane
species               Optional<String>            --species                              Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can speed up initial database connection as the registry does not have to load all available database aliases on the server. Default = "homo_sapiens"
assembly              Optional<String>            --assembly                             (-a) Select the assembly version to use if more than one available. If using the cache, you must
                                                                                                         have the appropriate assembly's cache file installed. If not specified and you have only 1 assembly
                                                                                                         version installed, this will be chosen by default. Default = use found assembly version
inputData             Optional<String>            --input_data                           (--id) Raw input data as a string. May be used, for example, to input a single rsID or HGVS notation quickly to vep: --input_data rs699
format                Optional<String>            --format                               Input file format - one of "ensembl", "vcf", "hgvs", "id", "region", "spdi". By default, VEP auto-detects the input file format. Using this option you can specify the input file is Ensembl, VCF, IDs, HGVS, SPDI or region format. Can use compressed version (gzipped) of any file format listed above. Auto-detects format by default
forceOverwrite        Optional<Boolean>           --force_overwrite                      (--force) By default, VEP will fail with an error if the output file already exists. You can force the overwrite of the existing file by using this flag. Not used by default
statsFile             Optional<String>            --stats_file                           (--sf) Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"
noStats               Optional<Boolean>           --no_stats                             Don't generate a stats file. Provides marginal gains in run time.
statsText             Optional<Boolean>           --stats_text                           Generate a plain text stats file in place of the HTML.
warningFile           Optional<Filename>          --warning_file                         File name to write warnings and errors to. Default = STDERR (standard error)
maxSvSize             Optional<Boolean>           --max_sv_size                          Extend the maximum Structural Variant size VEP can process.
noCheckVariantsOrder  Optional<Boolean>           --no_check_variants_order              Permit the use of unsorted input files. However running VEP on unsorted input files slows down the tool and requires more memory.
fork                  Optional<Integer>           --fork                                 Enable forking, using the specified number of forks. Forking can dramatically improve runtime. Not used by default
custom                Optional<Array<BedTABIX>>   --custom                               Add custom annotation to the output. Files must be tabix indexed or in the bigWig format. Multiple files can be specified by supplying the --custom flag multiple times. See https://asia.ensembl.org/info/docs/tools/vep/script/vep_custom.html for full details. Not used by default
gff                   Optional<File>              --gff                                  Use GFF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
gtf                   Optional<File>              --gtf                                  Use GTF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
bam                   Optional<BAM>               --bam                                  ADVANCED Use BAM file of sequence alignments to correct transcript models not derived from reference genome sequence. Used to correct RefSeq transcript models. Enables --use_transcript_ref; add --use_given_ref to override this behaviour. Not used by default
useTranscriptRef      Optional<Boolean>           --use_transcript_ref                   By default VEP uses the reference allele provided in the input file to calculate consequences for the provided alternate allele(s). Use this flag to force VEP to replace the provided reference allele with sequence derived from the overlapped transcript. This is especially relevant when using the RefSeq cache, see documentation for more details. The GIVEN_REF and USED_REF fields are set in the output to indicate any change. Not used by default
useGivenRef           Optional<Boolean>           --use_given_ref                        Using --bam or a BAM-edited RefSeq cache by default enables --use_transcript_ref; add this flag to override this behaviour and use the provided reference allele from the input. Not used by default
customMultiAllelic    Optional<Boolean>           --custom_multi_allelic                 By default, comma separated lists found within the INFO field of custom annotation VCFs are assumed to be allele specific. For example, a variant with allele_string A/G/C with associated custom annotation "single,double,triple" will associate triple with C, double with G and single with A. This flag instructs VEP to return all annotations for all alleles. Not used by default
tab                   Optional<Boolean>           --tab                                  Writes output in tab-delimited format. Not used by default
json                  Optional<Boolean>           --json                                 Writes output in JSON format. Not used by default
compressOutput        Optional<String>            --compress_output                      Writes output compressed using either gzip or bgzip. Not used by default
fields                Optional<Array<String>>     --fields                               Configure the output format using a comma separated list of fields.
                                                                                         Can only be used with tab (--tab) or VCF format (--vcf) output.
                                                                                         For the tab format output, the selected fields may be those present in the default output columns, or
                                                                                         any of those that appear in the Extra column (including those added by plugins or custom annotations).
                                                                                         Output remains tab-delimited. For the VCF format output, the selected fields are those present within the ""CSQ"" INFO field.

                                                                                         Example of command for the tab output:

                                                                                         --tab --fields ""Uploaded_variation,Location,Allele,Gene""
                                                                                         Example of command for the VCF format output:

                                                                                         --vcf --fields ""Allele,Consequence,Feature_type,Feature""
                                                                                         Not used by default
minimal               Optional<Boolean>           --minimal                              Convert alleles to their most minimal representation before consequence calculation i.e. sequence that is identical between each pair of reference and alternate alleles is trimmed off from both ends, with coordinates adjusted accordingly. Note this may lead to discrepancies between input coordinates and coordinates reported by VEP relative to transcript sequences; to avoid issues, use --allele_number and/or ensure that your input variants have unique identifiers. The MINIMISED flag is set in the VEP output where relevant. Not used by default
variantClass          Optional<Boolean>           --variant_class                        Output the Sequence Ontology variant class. Not used by default
sift                  Optional<String>            --sift                                 Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both. Not used by default
polyphen              Optional<String>            --polyphen                             Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both. VEP uses the humVar score by default - use --humdiv to retrieve the humDiv score. Not used by default
humdiv                Optional<Boolean>           --humdiv                               Human only Retrieve the humDiv PolyPhen prediction instead of the default humVar. Not used by default
nearest               Optional<String>            --nearest                              Retrieve the transcript or gene with the nearest protein-coding transcription start site
                                                                                                         (TSS) to each input variant. Use ""transcript"" to retrieve the transcript stable ID, ""gene"" to
                                                                                                         retrieve the gene stable ID, or ""symbol"" to retrieve the gene symbol. Note that the nearest
                                                                                                         TSS may not belong to a transcript that overlaps the input variant, and more than one may be
                                                                                                         reported in the case where two are equidistant from the input coordinates.

                                                                                                     Currently only available when using a cache annotation source, and requires the Set::IntervalTree perl module.
                                                                                                     Not used by default
distance              Optional<Array<Integer>>    --distance                             Modify the distance up and/or downstream between a variant and a transcript for which VEP will assign the upstream_gene_variant or downstream_gene_variant consequences. Giving one distance will modify both up- and downstream distances; prodiving two separated by commas will set the up- (5') and down - (3') stream distances respectively. Default: 5000
overlaps              Optional<Boolean>           --overlaps                             Report the proportion and length of a transcript overlapped by a structural variant in VCF format.
genePhenotype         Optional<Boolean>           --gene_phenotype                       Indicates if the overlapped gene is associated with a phenotype, disease or trait. See list of phenotype sources. Not used by default
regulatory            Optional<Boolean>           --regulatory                           Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature. Not used by default
cellType              Optional<Boolean>           --cell_type                            Report only regulatory regions that are found in the given cell type(s). Can be a single cell type or a comma-separated list. The functional type in each cell type is reported under CELL_TYPE in the output. To retrieve a list of cell types, use --cell_type list. Not used by default
individual            Optional<Array<String>>     --individual                           Consider only alternate alleles present in the genotypes of the specified individual(s). May be a single individual, a comma-separated list or "all" to assess all individuals separately. Individual variant combinations homozygous for the given reference allele will not be reported. Each individual and variant combination is given on a separate line of output. Only works with VCF files containing individual genotype data; individual IDs are taken from column headers. Not used by default
phased                Optional<Boolean>           --phased                               Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data. Not used by default
alleleNumber          Optional<Boolean>           --allele_number                        Identify allele number from VCF input, where 1 = first ALT allele, 2 = second ALT allele etc. Useful when using --minimal Not used by default
showRefAllele         Optional<Boolean>           --show_ref_allele                      Adds the reference allele in the output. Mainly useful for the VEP "default" and tab-delimited output formats. Not used by default
totalLength           Optional<Boolean>           --total_length                         Give cDNA, CDS and protein positions as Position/Length. Not used by default
numbers               Optional<Boolean>           --numbers                              Adds affected exon and intron numbering to to output. Format is Number/Total. Not used by default
noEscape              Optional<Boolean>           --no_escape                            Don't URI escape HGVS strings. Default = escape
keepCsq               Optional<Boolean>           --keep_csq                             Don't overwrite existing CSQ entry in VCF INFO field. Overwrites by default
vcfInfoField          Optional<String>            --vcf_info_field                       Change the name of the INFO key that VEP write the consequences to in its VCF output. Use "ANN" for compatibility with other tools such as snpEff. Default: CSQ
terms                 Optional<String>            --terms                                (-t) The type of consequence terms to output. The Ensembl terms are described here. The Sequence Ontology is a joint effort by genome annotation centres to standardise descriptions of biological sequences. Default = "SO"
noHeaders             Optional<Boolean>           --no_headers                           Don't write header lines in output files. Default = add headers
hgvs                  Optional<Boolean>           --hgvs                                 Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. HGVS notations given on Ensembl identifiers are versioned. Not used by default
hgvsg                 Optional<Boolean>           --hgvsg                                Add genomic HGVS nomenclature based on the input chromosome name. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. Not used by default
shiftHgvs             Optional<Boolean>           --shift_hgvs                           Enable or disable 3' shifting of HGVS notations. When enabled, this causes ambiguous insertions or deletions (typically in repetetive sequence tracts) to be "shifted" to their most 3' possible coordinates (relative to the transcript sequence and strand) before the HGVS notations are calculated; the flag HGVS_OFFSET is set to the number of bases by which the variant has shifted, relative to the input genomic coordinates. Disabling retains the original input coordinates of the variant. Default: 1 (shift)
transcriptVersion     Optional<Boolean>           --transcript_version                   Add version numbers to Ensembl transcript identifiers
protein               Optional<Boolean>           --protein                              Add the Ensembl protein identifier to the output where appropriate. Not used by default
symbol                Optional<Boolean>           --symbol                               Adds the gene symbol (e.g. HGNC) (where available) to the output. Not used by default
ccds                  Optional<Boolean>           --ccds                                 Adds the CCDS transcript identifer (where available) to the output. Not used by default
uniprot               Optional<Boolean>           --uniprot                              Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. Not used by default
tsl                   Optional<Boolean>           --tsl                                  Adds the transcript support level for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
appris                Optional<Boolean>           --appris                               Adds the APPRIS isoform annotation for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
canonical             Optional<Boolean>           --canonical                            Adds a flag indicating if the transcript is the canonical transcript for the gene. Not used by default
mane                  Optional<Boolean>           --mane                                 Adds a flag indicating if the transcript is the MANE Select transcript for the gene. Not used by default. Note: Only available for human on the GRCh38 assembly
biotype               Optional<Boolean>           --biotype                              Adds the biotype of the transcript or regulatory feature. Not used by default
domains               Optional<Boolean>           --domains                              Adds names of overlapping protein domains to output. Not used by default
xrefRefseq            Optional<Boolean>           --xref_refseq                          Output aligned RefSeq mRNA identifier for transcript. Not used by default. Note: The RefSeq and Ensembl transcripts aligned in this way MAY NOT, AND FREQUENTLY WILL NOT, match exactly in sequence, exon structure and protein product
synonyms              Optional<tsv>               --synonyms                             Load a file of chromosome synonyms. File should be tab-delimited with the primary identifier in column 1 and the synonym in column 2. Synonyms allow different chromosome identifiers to be used in the input file and any annotation source (cache, database, GFF, custom file, FASTA file). Not used by default
checkExisting         Optional<Boolean>           --check_existing                       Checks for the existence of known variants that are co-located with your input. By default the alleles are compared and variants on an allele-specific basis - to compare only coordinates, use --no_check_alleles.

                                                                                                     Some databases may contain variants with unknown (null) alleles and these are included by default; to exclude them use --exclude_null_alleles.

                                                                                                     See this page for more details.

                                                                                                     Not used by default
checkSvs              Optional<Boolean>           --check_svs                            Checks for the existence of structural variants that overlap your input. Currently requires database access. Not used by default
clinSigAllele         Optional<Boolean>           --clin_sig_allele                      Return allele specific clinical significance. Setting this option to 0 will provide all known clinical significance values at the given locus. Default: 1 (Provide allele-specific annotations)
excludeNullAlleles    Optional<Boolean>           --exclude_null_alleles                 Do not include variants with unknown alleles when checking for co-located variants. Our human database contains variants from HGMD and COSMIC for which the alleles are not publically available; by default these are included when using --check_existing, use this flag to exclude them. Not used by default
noCheckAlleles        Optional<Boolean>           --no_check_alleles                     When checking for existing variants, by default VEP only reports a co-located variant if none of the input alleles are novel. For example, if your input variant has alleles A/G, and an existing co-located variant has alleles A/C, the co-located variant will not be reported.

                                                                                                     Strand is also taken into account - in the same example, if the input variant has alleles T/G but on the negative strand, then the co-located variant will be reported since its alleles match the reverse complement of input variant.

                                                                                                     Use this flag to disable this behaviour and compare using coordinates alone. Not used by default
af                    Optional<Boolean>           --af                                   Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output. For this and all --af_* flags, the frequency reported is for the input allele only, not necessarily the non-reference or derived allele. Not used by default
maxAf                 Optional<Boolean>           --max_af                               Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD. Not used by default
af1kg                 Optional<String>            --af_1kg                               Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output. Must be used with --cache. Not used by default
afEsp                 Optional<Boolean>           --af_esp                               Include allele frequency from NHLBI-ESP populations. Must be used with --cache. Not used by default
afGnomad              Optional<Boolean>           --af_gnomad                            Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included; to retrieve data from the additional genomes data set, see this guide. Must be used with --cache Not used by default
afExac                Optional<Boolean>           --af_exac                              Include allele frequency from ExAC project populations. Must be used with --cache. Not used by default. Note: ExAC data has been superceded by gnomAD. This flag remains for those wishing to use older cache versions containing ExAC data.
pubmed                Optional<Boolean>           --pubmed                               Report Pubmed IDs for publications that cite existing variant. Must be used with --cache. Not used by default
failed                Optional<Boolean>           --failed                               When checking for co-located variants, by default VEP will exclude variants that have been flagged as failed. Set this flag to include such variants. Default: 0 (exclude)
gencodeBasic          Optional<Boolean>           --gencode_basic                        Limit your analysis to transcripts belonging to the GENCODE basic set. This set has fragmented or problematic transcripts removed. Not used by default
excludePredicted      Optional<Boolean>           --exclude_predicted                    When using the RefSeq or merged cache, exclude predicted transcripts (i.e. those with identifiers beginning with "XM_" or "XR_").
transcriptFilter      Optional<Boolean>           --transcript_filter                    ADVANCED Filter transcripts according to any arbitrary set of rules. Uses similar notation to filter_vep.

                                                                                                     You may filter on any key defined in the root of the transcript object; most commonly this will be ""stable_id"":

                                                                                                     --transcript_filter ""stable_id match N[MR]_""
checkRef              Optional<Boolean>           --check_ref                            Force VEP to check the supplied reference allele against the sequence stored in the Ensembl Core database or supplied FASTA file. Lines that do not match are skipped. Not used by default
lookupRef             Optional<Boolean>           --lookup_ref                           Force overwrite the supplied reference allele with the sequence stored in the Ensembl Core database or supplied FASTA file. Not used by default
dontSkip              Optional<Boolean>           --dont_skip                            Don't skip input variants that fail validation, e.g. those that fall on unrecognised sequences. Combining --check_ref with --dont_skip will add a CHECK_REF output field when the given reference does not match the underlying reference sequence.
allowNonVariant       Optional<Boolean>           --allow_non_variant                    When using VCF format as input and output, by default VEP will skip non-variant lines of input (where the ALT allele is null). Enabling this option the lines will be printed in the VCF output with no consequence data added.
chr                   Optional<Array<String>>     --chr                                  Select a subset of chromosomes to analyse from your file. Any data not on this chromosome in the input will be skipped. The list can be comma separated, with "-" characters representing an interval. For example, to include chromosomes 1, 2, 3, 10 and X you could use --chr 1-3,10,X Not used by default
codingOnly            Optional<Boolean>           --coding_only                          Only return consequences that fall in the coding regions of transcripts. Not used by default
noIntergenic          Optional<Boolean>           --no_intergenic                        Do not include intergenic consequences in the output. Not used by default
pick                  Optional<Boolean>           --pick                                 Pick once line or block of consequence data per variant, including transcript-specific columns. Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with --pick_order. This is the best method to use if you are interested only in one consequence per variant. Not used by default
pickAllele            Optional<Boolean>           --pick_allele                          Like --pick, but chooses one line or block of consequence data per variant allele. Will only differ in behaviour from --pick when the input variant has multiple alternate alleles. Not used by default
perGene               Optional<Boolean>           --per_gene                             Output only the most severe consequence per gene. The transcript selected is arbitrary if more than one has the same predicted consequence. Uses the same ranking system as --pick. Not used by default
pickAlleleGene        Optional<Boolean>           --pick_allele_gene                     Like --pick_allele, but chooses one line or block of consequence data per variant allele and gene combination. Not used by default
flagPick              Optional<Boolean>           --flag_pick                            As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
flagPickAllele        Optional<Boolean>           --flag_pick_allele                     As per --pick_allele, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
flagPickAlleleGene    Optional<Boolean>           --flag_pick_allele_gene                As per --pick_allele_gene, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
pickOrder             Optional<Array<String>>     --pick_order                           Customise the order of criteria (and the list of criteria) applied when choosing a block of annotation data with one of the following options: --pick, --pick_allele, --per_gene, --pick_allele_gene, --flag_pick, --flag_pick_allele, --flag_pick_allele_gene. See this page for the default order.
                                                                                                     Valid criteria are: [ canonical appris tsl biotype ccds rank length mane ]. e.g.:

                                                                                                     --pick --pick_order tsl,appris,rank
mostSevere            Optional<Boolean>           --most_severe                          Output only the most severe consequence per variant. Transcript-specific columns will be left blank. Consequence ranks are given in this table. To include regulatory consequences, use the --regulatory option in combination with this flag. Not used by default
summary               Optional<Boolean>           --summary                              Output only a comma-separated list of all observed consequences per variant. Transcript-specific columns will be left blank. Not used by default
filterCommon          Optional<Boolean>           --filter_common                        Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters. Not used by default
checkFrequency        Optional<Boolean>           --check_frequency                      Turns on frequency filtering. Use this to include or exclude variants based on the frequency of co-located existing variants in the Ensembl Variation database. You must also specify all of the --freq_* flags below. Frequencies used in filtering are added to the output under the FREQS key in the Extra field. Not used by default
freqPop               Optional<String>            --freq_pop                             Name of the population to use in frequency filter. This must be one of the following: (1KG_ALL, 1KG_AFR, 1KG_AMR, 1KG_EAS, 1KG_EUR, 1KG_SAS, AA, EA, gnomAD, gnomAD_AFR, gnomAD_AMR, gnomAD_ASJ, gnomAD_EAS, gnomAD_FIN, gnomAD_NFE, gnomAD_OTH, gnomAD_SAS)
freqFreq              Optional<Float>             --freq_freq                            Allele frequency to use for filtering. Must be a float value between 0 and 1
freqGtLt              Optional<String>            --freq_gt_lt                           Specify whether the frequency of the co-located variant must be greater than (gt) or less than (lt) the value specified with --freq_freq
freqFilter            Optional<String>            --freq_filter                          Specify whether to exclude or include only variants that pass the frequency filter
cache                 Optional<Boolean>           --cache                                Enables use of the cache. Add --refseq or --merged to use the refseq or merged cache.
cacheDir              Optional<Directory>         --dir                                  Specify the base cache/plugin directory to use. Default = "$HOME/.vep/"
dirCache              Optional<Directory>         --dir_cache                            Specify the cache directory to use. Default = "$HOME/.vep/"
dirPlugins            Optional<Directory>         --dir_plugins                          Specify the plugin directory to use. Default = "$HOME/.vep/"
offline               Optional<Boolean>           --offline                              Enable offline mode. No database connections will be made, and a cache file or GFF/GTF file is required for annotation. Add --refseq to use the refseq cache (if installed). Not used by default
fasta                 Optional<FastaWithIndexes>  --fasta                                (--fa) Specify a FASTA file or a directory containing FASTA files to use to look up reference sequence. The first time you run VEP with this parameter an index will be built which can take a few minutes. This is required if fetching HGVS annotations (--hgvs) or checking reference sequences (--check_ref) in offline mode (--offline), and optional with some performance increase in cache mode (--cache). See documentation for more details. Not used by default
refseq                Optional<Boolean>           --refseq                               Specify this option if you have installed the RefSeq cache in order for VEP to pick up the alternate cache directory. This cache contains transcript objects corresponding to RefSeq transcripts. Consequence output will be given relative to these transcripts in place of the default Ensembl transcripts (see documentation)
merged                Optional<Boolean>           --merged                               Use the merged Ensembl and RefSeq cache. Consequences are flagged with the SOURCE of each transcript used.
cacheVersion          Optional<Boolean>           --cache_version                        Use a different cache version than the assumed default (the VEP version). This should be used with Ensembl Genomes caches since their version numbers do not match Ensembl versions. For example, the VEP/Ensembl version may be 88 and the Ensembl Genomes version 35. Not used by default
showCacheInfo         Optional<Boolean>           --show_cache_info                      Show source version information for selected cache and quit
bufferSize            Optional<Integer>           --buffer_size                          Sets the internal buffer size, corresponding to the number of variants that are read in to memory simultaneously. Set this lower to use less memory at the expense of longer run time, and higher to use more memory with a faster run time. Default = 5000
====================  ==========================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
