:orphan:

Vep (Database)
=============================

``vep_database`` · *1 contributor · 1 version*

No documentation was provided: `contribute one <https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics>`_


Quickstart
-----------

    .. code-block:: python

       from janis_bioinformatics.tools.ensembl.vep.v98_3 import VepDatabase_98_3

       wf = WorkflowBuilder("myworkflow")

       wf.step(
           "vep_database_step",
           VepDatabase_98_3(
               inputFile=None,
               vcf=None,
           )
       )
       wf.output("out", source=vep_database_step.out)
       wf.output("out_stdout", source=vep_database_step.out_stdout)
       wf.output("out_stats", source=vep_database_step.out_stats)
       wf.output("out_warnings", source=vep_database_step.out_warnings)
    

*OR*

1. `Install Janis </tutorials/tutorial0.html>`_

2. Ensure Janis is configured to work with Docker or Singularity.

3. Ensure all reference files are available:

.. note:: 

   More information about these inputs are available `below <#additional-configuration-inputs>`_.



4. Generate user input files for vep_database:

.. code-block:: bash

   # user inputs
   janis inputs vep_database > inputs.yaml



**inputs.yaml**

.. code-block:: yaml

       inputFile: inputFile.vcf.gz




5. Run vep_database with:

.. code-block:: bash

   janis run [...run options] \
       --inputs inputs.yaml \
       vep_database





Information
------------

:ID: ``vep_database``
:URL: *No URL to the documentation was provided*
:Versions: 98.3
:Container: quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0
:Authors: Michael Franklin
:Citations: None
:Created: 2020-02-25
:Updated: 2020-05-07


Outputs
-----------

============  ======================  ===============
name          type                    documentation
============  ======================  ===============
out           Optional<Gzipped<VCF>>
out_stdout    stdout<TextFile>
out_stats     Optional<File>
out_warnings  Optional<File>
============  ======================  ===============


Additional configuration (inputs)
---------------------------------

====================  =============================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
name                  type                           prefix                     position    documentation
====================  =============================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================
inputFile             Gzipped<VCF>                   --input_file                           Input file name. Can use compressed file (gzipped).
vcf                   Boolean                        --vcf                                  Writes output in VCF format. Consequences are added in the INFO field of the VCF file, using the key "CSQ". Data fields are encoded separated by "|"; the order of fields is written in the VCF header. Output fields in the "CSQ" INFO field can be selected by using --fields. If the input format was VCF, the file will remain unchanged save for the addition of the CSQ field (unless using any filtering). Custom data added with --custom are added as separate fields, using the key specified for each data file. Commas in fields are replaced with ampersands (&) to preserve VCF format.
outputFilename        Optional<Filename>             --output_file                          (-o) Output file name. Results can write to STDOUT by specifying  as the output file name - this will force quiet mode. Default = "variant_effect_output.txt"
help                  Optional<Boolean>              --help                                 Display help message and quit
quiet                 Optional<Boolean>              --quiet                                (-q) Suppress warning messages.Not used by default
verbose               Optional<Boolean>              --verbose                              (-v) Print out a bit more information while running. Not used by default
config                Optional<File>                 --config                               Load configuration options from a config file. The config file should consist of whitespace-separated pairs of option names and settings e.g.:

                                                                                                        output_file   my_output.txt
                                                                                                        species       mus_musculus
                                                                                                        format        vcf
                                                                                                        host          useastdb.ensembl.org

                                                                                                        A config file can also be implicitly read; save the file as $HOME/.vep/vep.ini (or equivalent directory if
                                                                                                        using --dir). Any options in this file will be overridden by those specified in a config file using --config,
                                                                                                        and in turn by any options specified on the command line. You can create a quick version file of this by
                                                                                                        setting the flags as normal and running VEP in verbose (-v) mode. This will output lines that can be copied
                                                                                                        to a config file that can be loaded in on the next run using --config. Not used by default
everything            Optional<Boolean>              --everything                           (-e) Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --variant_class, --mane
species               Optional<String>               --species                              Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can speed up initial database connection as the registry does not have to load all available database aliases on the server. Default = "homo_sapiens"
assembly              Optional<String>               --assembly                             (-a) Select the assembly version to use if more than one available. If using the cache, you must
                                                                                                            have the appropriate assembly's cache file installed. If not specified and you have only 1 assembly
                                                                                                            version installed, this will be chosen by default. Default = use found assembly version
inputData             Optional<String>               --input_data                           (--id) Raw input data as a string. May be used, for example, to input a single rsID or HGVS notation quickly to vep: --input_data rs699
format                Optional<String>               --format                               Input file format - one of "ensembl", "vcf", "hgvs", "id", "region", "spdi". By default, VEP auto-detects the input file format. Using this option you can specify the input file is Ensembl, VCF, IDs, HGVS, SPDI or region format. Can use compressed version (gzipped) of any file format listed above. Auto-detects format by default
forceOverwrite        Optional<Boolean>              --force_overwrite                      (--force) By default, VEP will fail with an error if the output file already exists. You can force the overwrite of the existing file by using this flag. Not used by default
statsFile             Optional<String>               --stats_file                           (--sf) Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"
noStats               Optional<Boolean>              --no_stats                             Don't generate a stats file. Provides marginal gains in run time.
statsText             Optional<Boolean>              --stats_text                           Generate a plain text stats file in place of the HTML.
warningFile           Optional<Filename>             --warning_file                         File name to write warnings and errors to. Default = STDERR (standard error)
maxSvSize             Optional<Boolean>              --max_sv_size                          Extend the maximum Structural Variant size VEP can process.
noCheckVariantsOrder  Optional<Boolean>              --no_check_variants_order              Permit the use of unsorted input files. However running VEP on unsorted input files slows down the tool and requires more memory.
fork                  Optional<Integer>              --fork                                 Enable forking, using the specified number of forks. Forking can dramatically improve runtime. Not used by default
custom                Optional<Array<Gzipped<bed>>>  --custom                               Add custom annotation to the output. Files must be tabix indexed or in the bigWig format. Multiple files can be specified by supplying the --custom flag multiple times. See https://asia.ensembl.org/info/docs/tools/vep/script/vep_custom.html for full details. Not used by default
gff                   Optional<File>                 --gff                                  Use GFF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
gtf                   Optional<File>                 --gtf                                  Use GTF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
bam                   Optional<BAM>                  --bam                                  ADVANCED Use BAM file of sequence alignments to correct transcript models not derived from reference genome sequence. Used to correct RefSeq transcript models. Enables --use_transcript_ref; add --use_given_ref to override this behaviour. Not used by default
useTranscriptRef      Optional<Boolean>              --use_transcript_ref                   By default VEP uses the reference allele provided in the input file to calculate consequences for the provided alternate allele(s). Use this flag to force VEP to replace the provided reference allele with sequence derived from the overlapped transcript. This is especially relevant when using the RefSeq cache, see documentation for more details. The GIVEN_REF and USED_REF fields are set in the output to indicate any change. Not used by default
useGivenRef           Optional<Boolean>              --use_given_ref                        Using --bam or a BAM-edited RefSeq cache by default enables --use_transcript_ref; add this flag to override this behaviour and use the provided reference allele from the input. Not used by default
customMultiAllelic    Optional<Boolean>              --custom_multi_allelic                 By default, comma separated lists found within the INFO field of custom annotation VCFs are assumed to be allele specific. For example, a variant with allele_string A/G/C with associated custom annotation "single,double,triple" will associate triple with C, double with G and single with A. This flag instructs VEP to return all annotations for all alleles. Not used by default
tab                   Optional<Boolean>              --tab                                  Writes output in tab-delimited format. Not used by default
json                  Optional<Boolean>              --json                                 Writes output in JSON format. Not used by default
compressOutput        Optional<String>               --compress_output                      Writes output compressed using either gzip or bgzip. Not used by default
fields                Optional<Array<String>>        --fields                               Configure the output format using a comma separated list of fields.
                                                                                            Can only be used with tab (--tab) or VCF format (--vcf) output.
                                                                                            For the tab format output, the selected fields may be those present in the default output columns, or
                                                                                            any of those that appear in the Extra column (including those added by plugins or custom annotations).
                                                                                            Output remains tab-delimited. For the VCF format output, the selected fields are those present within the ""CSQ"" INFO field.

                                                                                            Example of command for the tab output:

                                                                                            --tab --fields ""Uploaded_variation,Location,Allele,Gene""
                                                                                            Example of command for the VCF format output:

                                                                                            --vcf --fields ""Allele,Consequence,Feature_type,Feature""
                                                                                            Not used by default
minimal               Optional<Boolean>              --minimal                              Convert alleles to their most minimal representation before consequence calculation i.e. sequence that is identical between each pair of reference and alternate alleles is trimmed off from both ends, with coordinates adjusted accordingly. Note this may lead to discrepancies between input coordinates and coordinates reported by VEP relative to transcript sequences; to avoid issues, use --allele_number and/or ensure that your input variants have unique identifiers. The MINIMISED flag is set in the VEP output where relevant. Not used by default
variantClass          Optional<Boolean>              --variant_class                        Output the Sequence Ontology variant class. Not used by default
sift                  Optional<String>               --sift                                 Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both. Not used by default
polyphen              Optional<String>               --polyphen                             Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both. VEP uses the humVar score by default - use --humdiv to retrieve the humDiv score. Not used by default
humdiv                Optional<Boolean>              --humdiv                               Human only Retrieve the humDiv PolyPhen prediction instead of the default humVar. Not used by default
nearest               Optional<String>               --nearest                              Retrieve the transcript or gene with the nearest protein-coding transcription start site
                                                                                                            (TSS) to each input variant. Use ""transcript"" to retrieve the transcript stable ID, ""gene"" to
                                                                                                            retrieve the gene stable ID, or ""symbol"" to retrieve the gene symbol. Note that the nearest
                                                                                                            TSS may not belong to a transcript that overlaps the input variant, and more than one may be
                                                                                                            reported in the case where two are equidistant from the input coordinates.

                                                                                                        Currently only available when using a cache annotation source, and requires the Set::IntervalTree perl module.
                                                                                                        Not used by default
distance              Optional<Array<Integer>>       --distance                             Modify the distance up and/or downstream between a variant and a transcript for which VEP will assign the upstream_gene_variant or downstream_gene_variant consequences. Giving one distance will modify both up- and downstream distances; prodiving two separated by commas will set the up- (5') and down - (3') stream distances respectively. Default: 5000
overlaps              Optional<Boolean>              --overlaps                             Report the proportion and length of a transcript overlapped by a structural variant in VCF format.
genePhenotype         Optional<Boolean>              --gene_phenotype                       Indicates if the overlapped gene is associated with a phenotype, disease or trait. See list of phenotype sources. Not used by default
regulatory            Optional<Boolean>              --regulatory                           Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature. Not used by default
cellType              Optional<Boolean>              --cell_type                            Report only regulatory regions that are found in the given cell type(s). Can be a single cell type or a comma-separated list. The functional type in each cell type is reported under CELL_TYPE in the output. To retrieve a list of cell types, use --cell_type list. Not used by default
individual            Optional<Array<String>>        --individual                           Consider only alternate alleles present in the genotypes of the specified individual(s). May be a single individual, a comma-separated list or "all" to assess all individuals separately. Individual variant combinations homozygous for the given reference allele will not be reported. Each individual and variant combination is given on a separate line of output. Only works with VCF files containing individual genotype data; individual IDs are taken from column headers. Not used by default
phased                Optional<Boolean>              --phased                               Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data. Not used by default
alleleNumber          Optional<Boolean>              --allele_number                        Identify allele number from VCF input, where 1 = first ALT allele, 2 = second ALT allele etc. Useful when using --minimal Not used by default
showRefAllele         Optional<Boolean>              --show_ref_allele                      Adds the reference allele in the output. Mainly useful for the VEP "default" and tab-delimited output formats. Not used by default
totalLength           Optional<Boolean>              --total_length                         Give cDNA, CDS and protein positions as Position/Length. Not used by default
numbers               Optional<Boolean>              --numbers                              Adds affected exon and intron numbering to to output. Format is Number/Total. Not used by default
noEscape              Optional<Boolean>              --no_escape                            Don't URI escape HGVS strings. Default = escape
keepCsq               Optional<Boolean>              --keep_csq                             Don't overwrite existing CSQ entry in VCF INFO field. Overwrites by default
vcfInfoField          Optional<String>               --vcf_info_field                       Change the name of the INFO key that VEP write the consequences to in its VCF output. Use "ANN" for compatibility with other tools such as snpEff. Default: CSQ
terms                 Optional<String>               --terms                                (-t) The type of consequence terms to output. The Ensembl terms are described here. The Sequence Ontology is a joint effort by genome annotation centres to standardise descriptions of biological sequences. Default = "SO"
noHeaders             Optional<Boolean>              --no_headers                           Don't write header lines in output files. Default = add headers
hgvs                  Optional<Boolean>              --hgvs                                 Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. HGVS notations given on Ensembl identifiers are versioned. Not used by default
hgvsg                 Optional<Boolean>              --hgvsg                                Add genomic HGVS nomenclature based on the input chromosome name. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. Not used by default
shiftHgvs             Optional<Boolean>              --shift_hgvs                           Enable or disable 3' shifting of HGVS notations. When enabled, this causes ambiguous insertions or deletions (typically in repetetive sequence tracts) to be "shifted" to their most 3' possible coordinates (relative to the transcript sequence and strand) before the HGVS notations are calculated; the flag HGVS_OFFSET is set to the number of bases by which the variant has shifted, relative to the input genomic coordinates. Disabling retains the original input coordinates of the variant. Default: 1 (shift)
transcriptVersion     Optional<Boolean>              --transcript_version                   Add version numbers to Ensembl transcript identifiers
protein               Optional<Boolean>              --protein                              Add the Ensembl protein identifier to the output where appropriate. Not used by default
symbol                Optional<Boolean>              --symbol                               Adds the gene symbol (e.g. HGNC) (where available) to the output. Not used by default
ccds                  Optional<Boolean>              --ccds                                 Adds the CCDS transcript identifer (where available) to the output. Not used by default
uniprot               Optional<Boolean>              --uniprot                              Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. Not used by default
tsl                   Optional<Boolean>              --tsl                                  Adds the transcript support level for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
appris                Optional<Boolean>              --appris                               Adds the APPRIS isoform annotation for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
canonical             Optional<Boolean>              --canonical                            Adds a flag indicating if the transcript is the canonical transcript for the gene. Not used by default
mane                  Optional<Boolean>              --mane                                 Adds a flag indicating if the transcript is the MANE Select transcript for the gene. Not used by default. Note: Only available for human on the GRCh38 assembly
biotype               Optional<Boolean>              --biotype                              Adds the biotype of the transcript or regulatory feature. Not used by default
domains               Optional<Boolean>              --domains                              Adds names of overlapping protein domains to output. Not used by default
xrefRefseq            Optional<Boolean>              --xref_refseq                          Output aligned RefSeq mRNA identifier for transcript. Not used by default. Note: The RefSeq and Ensembl transcripts aligned in this way MAY NOT, AND FREQUENTLY WILL NOT, match exactly in sequence, exon structure and protein product
synonyms              Optional<tsv>                  --synonyms                             Load a file of chromosome synonyms. File should be tab-delimited with the primary identifier in column 1 and the synonym in column 2. Synonyms allow different chromosome identifiers to be used in the input file and any annotation source (cache, database, GFF, custom file, FASTA file). Not used by default
checkExisting         Optional<Boolean>              --check_existing                       Checks for the existence of known variants that are co-located with your input. By default the alleles are compared and variants on an allele-specific basis - to compare only coordinates, use --no_check_alleles.

                                                                                                        Some databases may contain variants with unknown (null) alleles and these are included by default; to exclude them use --exclude_null_alleles.

                                                                                                        See this page for more details.

                                                                                                        Not used by default
checkSvs              Optional<Boolean>              --check_svs                            Checks for the existence of structural variants that overlap your input. Currently requires database access. Not used by default
clinSigAllele         Optional<Boolean>              --clin_sig_allele                      Return allele specific clinical significance. Setting this option to 0 will provide all known clinical significance values at the given locus. Default: 1 (Provide allele-specific annotations)
excludeNullAlleles    Optional<Boolean>              --exclude_null_alleles                 Do not include variants with unknown alleles when checking for co-located variants. Our human database contains variants from HGMD and COSMIC for which the alleles are not publically available; by default these are included when using --check_existing, use this flag to exclude them. Not used by default
noCheckAlleles        Optional<Boolean>              --no_check_alleles                     When checking for existing variants, by default VEP only reports a co-located variant if none of the input alleles are novel. For example, if your input variant has alleles A/G, and an existing co-located variant has alleles A/C, the co-located variant will not be reported.

                                                                                                        Strand is also taken into account - in the same example, if the input variant has alleles T/G but on the negative strand, then the co-located variant will be reported since its alleles match the reverse complement of input variant.

                                                                                                        Use this flag to disable this behaviour and compare using coordinates alone. Not used by default
af                    Optional<Boolean>              --af                                   Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output. For this and all --af_* flags, the frequency reported is for the input allele only, not necessarily the non-reference or derived allele. Not used by default
maxAf                 Optional<Boolean>              --max_af                               Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD. Not used by default
af1kg                 Optional<String>               --af_1kg                               Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output. Must be used with --cache. Not used by default
afEsp                 Optional<Boolean>              --af_esp                               Include allele frequency from NHLBI-ESP populations. Must be used with --cache. Not used by default
afGnomad              Optional<Boolean>              --af_gnomad                            Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included; to retrieve data from the additional genomes data set, see this guide. Must be used with --cache Not used by default
afExac                Optional<Boolean>              --af_exac                              Include allele frequency from ExAC project populations. Must be used with --cache. Not used by default. Note: ExAC data has been superceded by gnomAD. This flag remains for those wishing to use older cache versions containing ExAC data.
pubmed                Optional<Boolean>              --pubmed                               Report Pubmed IDs for publications that cite existing variant. Must be used with --cache. Not used by default
failed                Optional<Boolean>              --failed                               When checking for co-located variants, by default VEP will exclude variants that have been flagged as failed. Set this flag to include such variants. Default: 0 (exclude)
gencodeBasic          Optional<Boolean>              --gencode_basic                        Limit your analysis to transcripts belonging to the GENCODE basic set. This set has fragmented or problematic transcripts removed. Not used by default
excludePredicted      Optional<Boolean>              --exclude_predicted                    When using the RefSeq or merged cache, exclude predicted transcripts (i.e. those with identifiers beginning with "XM_" or "XR_").
transcriptFilter      Optional<Boolean>              --transcript_filter                    ADVANCED Filter transcripts according to any arbitrary set of rules. Uses similar notation to filter_vep.

                                                                                                        You may filter on any key defined in the root of the transcript object; most commonly this will be ""stable_id"":

                                                                                                        --transcript_filter ""stable_id match N[MR]_""
checkRef              Optional<Boolean>              --check_ref                            Force VEP to check the supplied reference allele against the sequence stored in the Ensembl Core database or supplied FASTA file. Lines that do not match are skipped. Not used by default
lookupRef             Optional<Boolean>              --lookup_ref                           Force overwrite the supplied reference allele with the sequence stored in the Ensembl Core database or supplied FASTA file. Not used by default
dontSkip              Optional<Boolean>              --dont_skip                            Don't skip input variants that fail validation, e.g. those that fall on unrecognised sequences. Combining --check_ref with --dont_skip will add a CHECK_REF output field when the given reference does not match the underlying reference sequence.
allowNonVariant       Optional<Boolean>              --allow_non_variant                    When using VCF format as input and output, by default VEP will skip non-variant lines of input (where the ALT allele is null). Enabling this option the lines will be printed in the VCF output with no consequence data added.
chr                   Optional<Array<String>>        --chr                                  Select a subset of chromosomes to analyse from your file. Any data not on this chromosome in the input will be skipped. The list can be comma separated, with "-" characters representing an interval. For example, to include chromosomes 1, 2, 3, 10 and X you could use --chr 1-3,10,X Not used by default
codingOnly            Optional<Boolean>              --coding_only                          Only return consequences that fall in the coding regions of transcripts. Not used by default
noIntergenic          Optional<Boolean>              --no_intergenic                        Do not include intergenic consequences in the output. Not used by default
pick                  Optional<Boolean>              --pick                                 Pick once line or block of consequence data per variant, including transcript-specific columns. Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with --pick_order. This is the best method to use if you are interested only in one consequence per variant. Not used by default
pickAllele            Optional<Boolean>              --pick_allele                          Like --pick, but chooses one line or block of consequence data per variant allele. Will only differ in behaviour from --pick when the input variant has multiple alternate alleles. Not used by default
perGene               Optional<Boolean>              --per_gene                             Output only the most severe consequence per gene. The transcript selected is arbitrary if more than one has the same predicted consequence. Uses the same ranking system as --pick. Not used by default
pickAlleleGene        Optional<Boolean>              --pick_allele_gene                     Like --pick_allele, but chooses one line or block of consequence data per variant allele and gene combination. Not used by default
flagPick              Optional<Boolean>              --flag_pick                            As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
flagPickAllele        Optional<Boolean>              --flag_pick_allele                     As per --pick_allele, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
flagPickAlleleGene    Optional<Boolean>              --flag_pick_allele_gene                As per --pick_allele_gene, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
pickOrder             Optional<Array<String>>        --pick_order                           Customise the order of criteria (and the list of criteria) applied when choosing a block of annotation data with one of the following options: --pick, --pick_allele, --per_gene, --pick_allele_gene, --flag_pick, --flag_pick_allele, --flag_pick_allele_gene. See this page for the default order.
                                                                                                        Valid criteria are: [ canonical appris tsl biotype ccds rank length mane ]. e.g.:

                                                                                                        --pick --pick_order tsl,appris,rank
mostSevere            Optional<Boolean>              --most_severe                          Output only the most severe consequence per variant. Transcript-specific columns will be left blank. Consequence ranks are given in this table. To include regulatory consequences, use the --regulatory option in combination with this flag. Not used by default
summary               Optional<Boolean>              --summary                              Output only a comma-separated list of all observed consequences per variant. Transcript-specific columns will be left blank. Not used by default
filterCommon          Optional<Boolean>              --filter_common                        Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters. Not used by default
checkFrequency        Optional<Boolean>              --check_frequency                      Turns on frequency filtering. Use this to include or exclude variants based on the frequency of co-located existing variants in the Ensembl Variation database. You must also specify all of the --freq_* flags below. Frequencies used in filtering are added to the output under the FREQS key in the Extra field. Not used by default
freqPop               Optional<String>               --freq_pop                             Name of the population to use in frequency filter. This must be one of the following: (1KG_ALL, 1KG_AFR, 1KG_AMR, 1KG_EAS, 1KG_EUR, 1KG_SAS, AA, EA, gnomAD, gnomAD_AFR, gnomAD_AMR, gnomAD_ASJ, gnomAD_EAS, gnomAD_FIN, gnomAD_NFE, gnomAD_OTH, gnomAD_SAS)
freqFreq              Optional<Float>                --freq_freq                            Allele frequency to use for filtering. Must be a float value between 0 and 1
freqGtLt              Optional<String>               --freq_gt_lt                           Specify whether the frequency of the co-located variant must be greater than (gt) or less than (lt) the value specified with --freq_freq
freqFilter            Optional<String>               --freq_filter                          Specify whether to exclude or include only variants that pass the frequency filter
caddReference         Optional<Array<Gzipped<VCF>>>
condelConfig          Optional<Directory>                                                   Directory containing CondelPlugin config, in format: '<dir>/condel_SP.conf'
dbnspReference        Optional<Gzipped<VCF>>
dbsnpColumns          Optional<Array<String>>
revelReference        Optional<Gzipped<VCF>>
custom1Reference      Optional<Gzipped<VCF>>
custom1Columns        Optional<Array<String>>
custom2Reference      Optional<Gzipped<VCF>>
custom2Columns        Optional<Array<String>>
database              Optional<Boolean>              --database                             Enable VEP to use local or remote databases.
host                  Optional<String>               --host                                 Manually define the database host to connect to. Users in the US may find connection and transfer speeds quicker using our East coast mirror, useastdb.ensembl.org. Default = "ensembldb.ensembl.org"
user                  Optional<String>               --user                                 (-u) Manually define the database username. Default = "anonymous"
password              Optional<String>               --password                             (--pass) Manually define the database password. Not used by default
port                  Optional<Integer>              --port                                 Manually define the database port. Default = 5306
genomes               Optional<Boolean>              --genomes                              Override the default connection settings with those for the Ensembl Genomes public MySQL server. Required when using any of the Ensembl Genomes species. Not used by default
isMultispecies        Optional<Boolean>              --is_multispecies                      Some of the Ensembl Genomes databases (mainly bacteria and protists) are composed of a collection of close species. It updates the database connection settings (i.e. the database name) if the value is set to 1. Default: 0
lrg                   Optional<Boolean>              --lrg                                  Map input variants to LRG coordinates (or to chromosome coordinates if given in LRG coordinates), and provide consequences on both LRG and chromosomal transcripts. Not used by default
dbVersion             Optional<String>               --db_version                           Force VEP to connect to a specific version of the Ensembl databases. Not recommended as there may be conflicts between software and database versions. Not used by default
registry              Optional<Filename>             --registry                             Defining a registry file overwrites other connection settings and uses those found in the specified registry file to connect. Not used by default
====================  =============================  =========================  ==========  =====================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================================

Workflow Description Language
------------------------------

.. code-block:: text

   version development

   task vep_database {
     input {
       Int? runtime_cpu
       Int? runtime_memory
       Int? runtime_seconds
       Int? runtime_disks
       File inputFile
       String? outputFilename
       Boolean? vcf
       Boolean? help
       Boolean? quiet
       Boolean? verbose
       File? config
       Boolean? everything
       String? species
       String? assembly
       String? inputData
       String? format
       Boolean? forceOverwrite
       String? statsFile
       Boolean? noStats
       Boolean? statsText
       String? warningFile
       Boolean? maxSvSize
       Boolean? noCheckVariantsOrder
       Int? fork
       Array[File]? custom
       Array[File]? custom_tbi
       File? gff
       File? gtf
       File? bam
       Boolean? useTranscriptRef
       Boolean? useGivenRef
       Boolean? customMultiAllelic
       Boolean? tab
       Boolean? json
       String? compressOutput
       Array[String]? fields
       Boolean? minimal
       Boolean? variantClass
       String? sift
       String? polyphen
       Boolean? humdiv
       String? nearest
       Array[Int]? distance
       Boolean? overlaps
       Boolean? genePhenotype
       Boolean? regulatory
       Boolean? cellType
       Array[String]? individual
       Boolean? phased
       Boolean? alleleNumber
       Boolean? showRefAllele
       Boolean? totalLength
       Boolean? numbers
       Boolean? noEscape
       Boolean? keepCsq
       String? vcfInfoField
       String? terms
       Boolean? noHeaders
       Boolean? hgvs
       Boolean? hgvsg
       Boolean? shiftHgvs
       Boolean? transcriptVersion
       Boolean? protein
       Boolean? symbol
       Boolean? ccds
       Boolean? uniprot
       Boolean? tsl
       Boolean? appris
       Boolean? canonical
       Boolean? mane
       Boolean? biotype
       Boolean? domains
       Boolean? xrefRefseq
       File? synonyms
       Boolean? checkExisting
       Boolean? checkSvs
       Boolean? clinSigAllele
       Boolean? excludeNullAlleles
       Boolean? noCheckAlleles
       Boolean? af
       Boolean? maxAf
       String? af1kg
       Boolean? afEsp
       Boolean? afGnomad
       Boolean? afExac
       Boolean? pubmed
       Boolean? failed
       Boolean? gencodeBasic
       Boolean? excludePredicted
       Boolean? transcriptFilter
       Boolean? checkRef
       Boolean? lookupRef
       Boolean? dontSkip
       Boolean? allowNonVariant
       Array[String]? chr
       Boolean? codingOnly
       Boolean? noIntergenic
       Boolean? pick
       Boolean? pickAllele
       Boolean? perGene
       Boolean? pickAlleleGene
       Boolean? flagPick
       Boolean? flagPickAllele
       Boolean? flagPickAlleleGene
       Array[String]? pickOrder
       Boolean? mostSevere
       Boolean? summary
       Boolean? filterCommon
       Boolean? checkFrequency
       String? freqPop
       Float? freqFreq
       String? freqGtLt
       String? freqFilter
       Array[File]? caddReference
       Array[File]? caddReference_tbi
       Directory? condelConfig
       File? dbnspReference
       File? dbnspReference_tbi
       Array[String]? dbsnpColumns
       File? revelReference
       File? revelReference_tbi
       File? custom1Reference
       File? custom1Reference_tbi
       Array[String]? custom1Columns
       File? custom2Reference
       File? custom2Reference_tbi
       Array[String]? custom2Columns
       Boolean? database
       String? host
       String? user
       String? password
       Int? port
       Boolean? genomes
       Boolean? isMultispecies
       Boolean? lrg
       String? dbVersion
       String? registry
     }
     command <<<
       set -e
       vep \
         --input_file '~{inputFile}' \
         --output_file '~{select_first([outputFilename, "~{basename(inputFile, ".vcf.gz")}.vcf"])}' \
         ~{if select_first([vcf, true]) then "--vcf" else ""} \
         ~{if (defined(help) && select_first([help])) then "--help" else ""} \
         ~{if (defined(quiet) && select_first([quiet])) then "--quiet" else ""} \
         ~{if (defined(verbose) && select_first([verbose])) then "--verbose" else ""} \
         ~{if defined(config) then ("--config '" + config + "'") else ""} \
         ~{if (defined(everything) && select_first([everything])) then "--everything" else ""} \
         ~{if defined(species) then ("--species '" + species + "'") else ""} \
         ~{if defined(assembly) then ("--assembly '" + assembly + "'") else ""} \
         ~{if defined(inputData) then ("--input_data '" + inputData + "'") else ""} \
         ~{if defined(format) then ("--format '" + format + "'") else ""} \
         ~{if (defined(forceOverwrite) && select_first([forceOverwrite])) then "--force_overwrite" else ""} \
         ~{if defined(select_first([statsFile, "variant_effect_output.txt_summary.html"])) then ("--stats_file '" + select_first([statsFile, "variant_effect_output.txt_summary.html"]) + "'") else ""} \
         ~{if (defined(noStats) && select_first([noStats])) then "--no_stats" else ""} \
         ~{if (defined(statsText) && select_first([statsText])) then "--stats_text" else ""} \
         --warning_file '~{select_first([warningFile, "generated-warning.txt"])}' \
         ~{if (defined(maxSvSize) && select_first([maxSvSize])) then "--max_sv_size" else ""} \
         ~{if (defined(noCheckVariantsOrder) && select_first([noCheckVariantsOrder])) then "--no_check_variants_order" else ""} \
         ~{if defined(select_first([fork, select_first([runtime_cpu, 1])])) then ("--fork " + select_first([fork, select_first([runtime_cpu, 1])])) else ''} \
         ~{if (defined(custom) && length(select_first([custom])) > 0) then "--custom '" + sep("' --custom '", select_first([custom])) + "'" else ""} \
         ~{if defined(gff) then ("--gff '" + gff + "'") else ""} \
         ~{if defined(gtf) then ("--gtf '" + gtf + "'") else ""} \
         ~{if defined(bam) then ("--bam '" + bam + "'") else ""} \
         ~{if (defined(useTranscriptRef) && select_first([useTranscriptRef])) then "--use_transcript_ref" else ""} \
         ~{if (defined(useGivenRef) && select_first([useGivenRef])) then "--use_given_ref" else ""} \
         ~{if (defined(customMultiAllelic) && select_first([customMultiAllelic])) then "--custom_multi_allelic" else ""} \
         ~{if (defined(tab) && select_first([tab])) then "--tab" else ""} \
         ~{if (defined(json) && select_first([json])) then "--json" else ""} \
         ~{if defined(select_first([compressOutput, "bgzip"])) then ("--compress_output '" + select_first([compressOutput, "bgzip"]) + "'") else ""} \
         ~{if (defined(fields) && length(select_first([fields])) > 0) then "--fields '" + sep("' '", select_first([fields])) + "'" else ""} \
         ~{if (defined(minimal) && select_first([minimal])) then "--minimal" else ""} \
         ~{if (defined(variantClass) && select_first([variantClass])) then "--variant_class" else ""} \
         ~{if defined(sift) then ("--sift '" + sift + "'") else ""} \
         ~{if defined(polyphen) then ("--polyphen '" + polyphen + "'") else ""} \
         ~{if (defined(humdiv) && select_first([humdiv])) then "--humdiv" else ""} \
         ~{if defined(nearest) then ("--nearest '" + nearest + "'") else ""} \
         ~{if (defined(distance) && length(select_first([distance])) > 0) then "--distance " + sep(",", select_first([distance])) else ""} \
         ~{if (defined(overlaps) && select_first([overlaps])) then "--overlaps" else ""} \
         ~{if (defined(genePhenotype) && select_first([genePhenotype])) then "--gene_phenotype" else ""} \
         ~{if (defined(regulatory) && select_first([regulatory])) then "--regulatory" else ""} \
         ~{if (defined(cellType) && select_first([cellType])) then "--cell_type" else ""} \
         ~{if (defined(individual) && length(select_first([individual])) > 0) then "--individual '" + sep("','", select_first([individual])) + "'" else ""} \
         ~{if (defined(phased) && select_first([phased])) then "--phased" else ""} \
         ~{if (defined(alleleNumber) && select_first([alleleNumber])) then "--allele_number" else ""} \
         ~{if (defined(showRefAllele) && select_first([showRefAllele])) then "--show_ref_allele" else ""} \
         ~{if (defined(totalLength) && select_first([totalLength])) then "--total_length" else ""} \
         ~{if (defined(numbers) && select_first([numbers])) then "--numbers" else ""} \
         ~{if (defined(noEscape) && select_first([noEscape])) then "--no_escape" else ""} \
         ~{if (defined(keepCsq) && select_first([keepCsq])) then "--keep_csq" else ""} \
         ~{if defined(vcfInfoField) then ("--vcf_info_field '" + vcfInfoField + "'") else ""} \
         ~{if defined(terms) then ("--terms '" + terms + "'") else ""} \
         ~{if (defined(noHeaders) && select_first([noHeaders])) then "--no_headers" else ""} \
         ~{if (defined(hgvs) && select_first([hgvs])) then "--hgvs" else ""} \
         ~{if (defined(hgvsg) && select_first([hgvsg])) then "--hgvsg" else ""} \
         ~{if (defined(shiftHgvs) && select_first([shiftHgvs])) then "--shift_hgvs" else ""} \
         ~{if (defined(transcriptVersion) && select_first([transcriptVersion])) then "--transcript_version" else ""} \
         ~{if (defined(protein) && select_first([protein])) then "--protein" else ""} \
         ~{if (defined(symbol) && select_first([symbol])) then "--symbol" else ""} \
         ~{if (defined(ccds) && select_first([ccds])) then "--ccds" else ""} \
         ~{if (defined(uniprot) && select_first([uniprot])) then "--uniprot" else ""} \
         ~{if (defined(tsl) && select_first([tsl])) then "--tsl" else ""} \
         ~{if (defined(appris) && select_first([appris])) then "--appris" else ""} \
         ~{if (defined(canonical) && select_first([canonical])) then "--canonical" else ""} \
         ~{if (defined(mane) && select_first([mane])) then "--mane" else ""} \
         ~{if (defined(biotype) && select_first([biotype])) then "--biotype" else ""} \
         ~{if (defined(domains) && select_first([domains])) then "--domains" else ""} \
         ~{if (defined(xrefRefseq) && select_first([xrefRefseq])) then "--xref_refseq" else ""} \
         ~{if defined(synonyms) then ("--synonyms '" + synonyms + "'") else ""} \
         ~{if (defined(checkExisting) && select_first([checkExisting])) then "--check_existing" else ""} \
         ~{if (defined(checkSvs) && select_first([checkSvs])) then "--check_svs" else ""} \
         ~{if (defined(clinSigAllele) && select_first([clinSigAllele])) then "--clin_sig_allele" else ""} \
         ~{if (defined(excludeNullAlleles) && select_first([excludeNullAlleles])) then "--exclude_null_alleles" else ""} \
         ~{if (defined(noCheckAlleles) && select_first([noCheckAlleles])) then "--no_check_alleles" else ""} \
         ~{if (defined(af) && select_first([af])) then "--af" else ""} \
         ~{if (defined(maxAf) && select_first([maxAf])) then "--max_af" else ""} \
         ~{if defined(af1kg) then ("--af_1kg '" + af1kg + "'") else ""} \
         ~{if (defined(afEsp) && select_first([afEsp])) then "--af_esp" else ""} \
         ~{if (defined(afGnomad) && select_first([afGnomad])) then "--af_gnomad" else ""} \
         ~{if (defined(afExac) && select_first([afExac])) then "--af_exac" else ""} \
         ~{if (defined(pubmed) && select_first([pubmed])) then "--pubmed" else ""} \
         ~{if (defined(failed) && select_first([failed])) then "--failed" else ""} \
         ~{if (defined(gencodeBasic) && select_first([gencodeBasic])) then "--gencode_basic" else ""} \
         ~{if (defined(excludePredicted) && select_first([excludePredicted])) then "--exclude_predicted" else ""} \
         ~{if (defined(transcriptFilter) && select_first([transcriptFilter])) then "--transcript_filter" else ""} \
         ~{if (defined(checkRef) && select_first([checkRef])) then "--check_ref" else ""} \
         ~{if (defined(lookupRef) && select_first([lookupRef])) then "--lookup_ref" else ""} \
         ~{if (defined(dontSkip) && select_first([dontSkip])) then "--dont_skip" else ""} \
         ~{if (defined(allowNonVariant) && select_first([allowNonVariant])) then "--allow_non_variant" else ""} \
         ~{if (defined(chr) && length(select_first([chr])) > 0) then "--chr '" + sep("','", select_first([chr])) + "'" else ""} \
         ~{if (defined(codingOnly) && select_first([codingOnly])) then "--coding_only" else ""} \
         ~{if (defined(noIntergenic) && select_first([noIntergenic])) then "--no_intergenic" else ""} \
         ~{if (defined(pick) && select_first([pick])) then "--pick" else ""} \
         ~{if (defined(pickAllele) && select_first([pickAllele])) then "--pick_allele" else ""} \
         ~{if (defined(perGene) && select_first([perGene])) then "--per_gene" else ""} \
         ~{if (defined(pickAlleleGene) && select_first([pickAlleleGene])) then "--pick_allele_gene" else ""} \
         ~{if (defined(flagPick) && select_first([flagPick])) then "--flag_pick" else ""} \
         ~{if (defined(flagPickAllele) && select_first([flagPickAllele])) then "--flag_pick_allele" else ""} \
         ~{if (defined(flagPickAlleleGene) && select_first([flagPickAlleleGene])) then "--flag_pick_allele_gene" else ""} \
         ~{if (defined(pickOrder) && length(select_first([pickOrder])) > 0) then "--pick_order '" + sep("','", select_first([pickOrder])) + "'" else ""} \
         ~{if (defined(mostSevere) && select_first([mostSevere])) then "--most_severe" else ""} \
         ~{if (defined(summary) && select_first([summary])) then "--summary" else ""} \
         ~{if (defined(filterCommon) && select_first([filterCommon])) then "--filter_common" else ""} \
         ~{if (defined(checkFrequency) && select_first([checkFrequency])) then "--check_frequency" else ""} \
         ~{if defined(freqPop) then ("--freq_pop '" + freqPop + "'") else ""} \
         ~{if defined(freqFreq) then ("--freq_freq " + freqFreq) else ''} \
         ~{if defined(freqGtLt) then ("--freq_gt_lt '" + freqGtLt + "'") else ""} \
         ~{if defined(freqFilter) then ("--freq_filter '" + freqFilter + "'") else ""} \
         ~{if (defined(database) && select_first([database])) then "--database" else ""} \
         ~{if defined(host) then ("--host '" + host + "'") else ""} \
         ~{if defined(user) then ("--user '" + user + "'") else ""} \
         ~{if defined(password) then ("--password '" + password + "'") else ""} \
         ~{if defined(port) then ("--port " + port) else ''} \
         ~{if (defined(genomes) && select_first([genomes])) then "--genomes" else ""} \
         ~{if (defined(isMultispecies) && select_first([isMultispecies])) then "--is_multispecies" else ""} \
         ~{if (defined(lrg) && select_first([lrg])) then "--lrg" else ""} \
         ~{if defined(dbVersion) then ("--db_version '" + dbVersion + "'") else ""} \
         --registry '~{select_first([registry, "generated"])}' \
         ~{if (defined(caddReference)) then ("--plugin CADD," + sep(",", select_first([caddReference]))) else ""} \
         ~{if (defined(condelConfig)) then "--plugin Condel,~{select_first([condelConfig])},b" else ""} \
         ~{if ((defined(dbnspReference) && defined(dbsnpColumns))) then "--plugin dbNSFP,~{select_first([dbnspReference])},~{sep(",", select_first([dbsnpColumns]))}" else ""} \
         ~{if (defined(revelReference)) then "--plugin REVEL,~{select_first([revelReference])}" else ""} \
         ~{if ((defined(custom1Reference) && defined(custom1Columns))) then "--custom ~{select_first([custom1Reference])},~{sep(",", select_first([custom1Columns]))}" else ""} \
         ~{if ((defined(custom2Reference) && defined(custom2Columns))) then "--custom ~{select_first([custom2Reference])},~{sep(",", select_first([custom2Columns]))}" else ""}
     >>>
     runtime {
       cpu: select_first([runtime_cpu, 1])
       disks: "local-disk ~{select_first([runtime_disks, 20])} SSD"
       docker: "quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0"
       duration: select_first([runtime_seconds, 86400])
       memory: "~{select_first([runtime_memory, 4])}G"
       preemptible: 2
     }
     output {
       File? out = select_first([outputFilename, "~{basename(inputFile, ".vcf.gz")}.vcf"])
       File? out_tbi = if defined(select_first([outputFilename, "~{basename(inputFile, ".vcf.gz")}.vcf"])) then (select_first([outputFilename, "~{basename(inputFile, ".vcf.gz")}.vcf"]) + ".tbi") else None
       File out_stdout = stdout()
       File? out_stats = select_first([statsFile, "variant_effect_output.txt_summary.html"])
       File? out_warnings = select_first([warningFile, "generated-warning.txt"])
     }
   }

Common Workflow Language
-------------------------

.. code-block:: text

   #!/usr/bin/env cwl-runner
   class: CommandLineTool
   cwlVersion: v1.2
   label: Vep (Database)
   doc: ''

   requirements:
   - class: ShellCommandRequirement
   - class: InlineJavascriptRequirement
   - class: DockerRequirement
     dockerPull: quay.io/biocontainers/ensembl-vep:98.3--pl526hecc5488_0

   inputs:
   - id: inputFile
     label: inputFile
     doc: Input file name. Can use compressed file (gzipped).
     type: File
     inputBinding:
       prefix: --input_file
   - id: outputFilename
     label: outputFilename
     doc: |-
       (-o) Output file name. Results can write to STDOUT by specifying  as the output file name - this will force quiet mode. Default = "variant_effect_output.txt"
     type:
     - string
     - 'null'
     default: generated.vcf
     inputBinding:
       prefix: --output_file
       valueFrom: $(inputs.inputFile.basename.replace(/.vcf.gz$/, "")).vcf
   - id: vcf
     label: vcf
     doc: |-
       Writes output in VCF format. Consequences are added in the INFO field of the VCF file, using the key "CSQ". Data fields are encoded separated by "|"; the order of fields is written in the VCF header. Output fields in the "CSQ" INFO field can be selected by using --fields. If the input format was VCF, the file will remain unchanged save for the addition of the CSQ field (unless using any filtering). Custom data added with --custom are added as separate fields, using the key specified for each data file. Commas in fields are replaced with ampersands (&) to preserve VCF format.
     type: boolean
     default: true
     inputBinding:
       prefix: --vcf
   - id: help
     label: help
     doc: Display help message and quit
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --help
   - id: quiet
     label: quiet
     doc: (-q) Suppress warning messages.Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --quiet
   - id: verbose
     label: verbose
     doc: (-v) Print out a bit more information while running. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --verbose
   - id: config
     label: config
     doc: |-
       Load configuration options from a config file. The config file should consist of whitespace-separated pairs of option names and settings e.g.:

                   output_file   my_output.txt
                   species       mus_musculus
                   format        vcf
                   host          useastdb.ensembl.org

                   A config file can also be implicitly read; save the file as $HOME/.vep/vep.ini (or equivalent directory if 
                   using --dir). Any options in this file will be overridden by those specified in a config file using --config, 
                   and in turn by any options specified on the command line. You can create a quick version file of this by 
                   setting the flags as normal and running VEP in verbose (-v) mode. This will output lines that can be copied 
                   to a config file that can be loaded in on the next run using --config. Not used by default
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --config
   - id: everything
     label: everything
     doc: |-
       (-e) Shortcut flag to switch on all of the following: --sift b, --polyphen b, --ccds, --uniprot, --hgvs, --symbol, --numbers, --domains, --regulatory, --canonical, --protein, --biotype, --uniprot, --tsl, --appris, --gene_phenotype --af, --af_1kg, --af_esp, --af_gnomad, --max_af, --pubmed, --variant_class, --mane
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --everything
   - id: species
     label: species
     doc: |-
       Species for your data. This can be the latin name e.g. "homo_sapiens" or any Ensembl alias e.g. "mouse". Specifying the latin name can speed up initial database connection as the registry does not have to load all available database aliases on the server. Default = "homo_sapiens"
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --species
   - id: assembly
     label: assembly
     doc: |-
       (-a) Select the assembly version to use if more than one available. If using the cache, you must 
                       have the appropriate assembly's cache file installed. If not specified and you have only 1 assembly 
                       version installed, this will be chosen by default. Default = use found assembly version
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --assembly
   - id: inputData
     label: inputData
     doc: |-
       (--id) Raw input data as a string. May be used, for example, to input a single rsID or HGVS notation quickly to vep: --input_data rs699
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --input_data
   - id: format
     label: format
     doc: |-
       Input file format - one of "ensembl", "vcf", "hgvs", "id", "region", "spdi". By default, VEP auto-detects the input file format. Using this option you can specify the input file is Ensembl, VCF, IDs, HGVS, SPDI or region format. Can use compressed version (gzipped) of any file format listed above. Auto-detects format by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --format
   - id: forceOverwrite
     label: forceOverwrite
     doc: |-
       (--force) By default, VEP will fail with an error if the output file already exists. You can force the overwrite of the existing file by using this flag. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --force_overwrite
   - id: statsFile
     label: statsFile
     doc: |-
       (--sf) Summary stats file name. This is an HTML file containing a summary of the VEP run - the file name must end ".htm" or ".html". Default = "variant_effect_output.txt_summary.html"
     type: string
     default: variant_effect_output.txt_summary.html
     inputBinding:
       prefix: --stats_file
   - id: noStats
     label: noStats
     doc: Don't generate a stats file. Provides marginal gains in run time.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_stats
   - id: statsText
     label: statsText
     doc: Generate a plain text stats file in place of the HTML.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --stats_text
   - id: warningFile
     label: warningFile
     doc: File name to write warnings and errors to. Default = STDERR (standard error)
     type:
     - string
     - 'null'
     default: generated-warning.txt
     inputBinding:
       prefix: --warning_file
   - id: maxSvSize
     label: maxSvSize
     doc: Extend the maximum Structural Variant size VEP can process.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --max_sv_size
   - id: noCheckVariantsOrder
     label: noCheckVariantsOrder
     doc: |-
       Permit the use of unsorted input files. However running VEP on unsorted input files slows down the tool and requires more memory.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_check_variants_order
   - id: fork
     label: fork
     doc: |-
       Enable forking, using the specified number of forks. Forking can dramatically improve runtime. Not used by default
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --fork
       valueFrom: $([inputs.runtime_cpu, 1].filter(function (inner) { return inner !=
         null })[0])
   - id: custom
     label: custom
     doc: |-
       Add custom annotation to the output. Files must be tabix indexed or in the bigWig format. Multiple files can be specified by supplying the --custom flag multiple times. See https://asia.ensembl.org/info/docs/tools/vep/script/vep_custom.html for full details. Not used by default
     type:
     - type: array
       inputBinding:
         prefix: --custom
       items: File
     - 'null'
     inputBinding: {}
   - id: gff
     label: gff
     doc: |-
       Use GFF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --gff
   - id: gtf
     label: gtf
     doc: |-
       Use GTF transcript annotations in [filename] as an annotation source. Requires a FASTA file of genomic sequence.Not used by default
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --gtf
   - id: bam
     label: bam
     doc: |-
       ADVANCED Use BAM file of sequence alignments to correct transcript models not derived from reference genome sequence. Used to correct RefSeq transcript models. Enables --use_transcript_ref; add --use_given_ref to override this behaviour. Not used by default
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --bam
   - id: useTranscriptRef
     label: useTranscriptRef
     doc: |-
       By default VEP uses the reference allele provided in the input file to calculate consequences for the provided alternate allele(s). Use this flag to force VEP to replace the provided reference allele with sequence derived from the overlapped transcript. This is especially relevant when using the RefSeq cache, see documentation for more details. The GIVEN_REF and USED_REF fields are set in the output to indicate any change. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_transcript_ref
   - id: useGivenRef
     label: useGivenRef
     doc: |-
       Using --bam or a BAM-edited RefSeq cache by default enables --use_transcript_ref; add this flag to override this behaviour and use the provided reference allele from the input. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --use_given_ref
   - id: customMultiAllelic
     label: customMultiAllelic
     doc: |-
       By default, comma separated lists found within the INFO field of custom annotation VCFs are assumed to be allele specific. For example, a variant with allele_string A/G/C with associated custom annotation "single,double,triple" will associate triple with C, double with G and single with A. This flag instructs VEP to return all annotations for all alleles. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --custom_multi_allelic
   - id: tab
     label: tab
     doc: Writes output in tab-delimited format. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --tab
   - id: json
     label: json
     doc: Writes output in JSON format. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --json
   - id: compressOutput
     label: compressOutput
     doc: Writes output compressed using either gzip or bgzip. Not used by default
     type: string
     default: bgzip
     inputBinding:
       prefix: --compress_output
   - id: fields
     label: fields
     doc: |-
       Configure the output format using a comma separated list of fields.
       Can only be used with tab (--tab) or VCF format (--vcf) output.
       For the tab format output, the selected fields may be those present in the default output columns, or 
       any of those that appear in the Extra column (including those added by plugins or custom annotations). 
       Output remains tab-delimited. For the VCF format output, the selected fields are those present within the ""CSQ"" INFO field.

       Example of command for the tab output:

       --tab --fields ""Uploaded_variation,Location,Allele,Gene""
       Example of command for the VCF format output:

       --vcf --fields ""Allele,Consequence,Feature_type,Feature""
       Not used by default
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --fields
   - id: minimal
     label: minimal
     doc: |-
       Convert alleles to their most minimal representation before consequence calculation i.e. sequence that is identical between each pair of reference and alternate alleles is trimmed off from both ends, with coordinates adjusted accordingly. Note this may lead to discrepancies between input coordinates and coordinates reported by VEP relative to transcript sequences; to avoid issues, use --allele_number and/or ensure that your input variants have unique identifiers. The MINIMISED flag is set in the VEP output where relevant. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --minimal
   - id: variantClass
     label: variantClass
     doc: Output the Sequence Ontology variant class. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --variant_class
   - id: sift
     label: sift
     doc: |-
       Species limited SIFT predicts whether an amino acid substitution affects protein function based on sequence homology and the physical properties of amino acids. VEP can output the prediction term, score or both. Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --sift
   - id: polyphen
     label: polyphen
     doc: |-
       Human only PolyPhen is a tool which predicts possible impact of an amino acid substitution on the structure and function of a human protein using straightforward physical and comparative considerations. VEP can output the prediction term, score or both. VEP uses the humVar score by default - use --humdiv to retrieve the humDiv score. Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --polyphen
   - id: humdiv
     label: humdiv
     doc: |-
       Human only Retrieve the humDiv PolyPhen prediction instead of the default humVar. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --humdiv
   - id: nearest
     label: nearest
     doc: |-
       Retrieve the transcript or gene with the nearest protein-coding transcription start site 
                       (TSS) to each input variant. Use ""transcript"" to retrieve the transcript stable ID, ""gene"" to 
                       retrieve the gene stable ID, or ""symbol"" to retrieve the gene symbol. Note that the nearest 
                       TSS may not belong to a transcript that overlaps the input variant, and more than one may be 
                       reported in the case where two are equidistant from the input coordinates.

                   Currently only available when using a cache annotation source, and requires the Set::IntervalTree perl module.
                   Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --nearest
   - id: distance
     label: distance
     doc: |-
       Modify the distance up and/or downstream between a variant and a transcript for which VEP will assign the upstream_gene_variant or downstream_gene_variant consequences. Giving one distance will modify both up- and downstream distances; prodiving two separated by commas will set the up- (5') and down - (3') stream distances respectively. Default: 5000
     type:
     - type: array
       items: int
     - 'null'
     inputBinding:
       prefix: --distance
       itemSeparator: ','
   - id: overlaps
     label: overlaps
     doc: |-
       Report the proportion and length of a transcript overlapped by a structural variant in VCF format.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --overlaps
   - id: genePhenotype
     label: genePhenotype
     doc: |-
       Indicates if the overlapped gene is associated with a phenotype, disease or trait. See list of phenotype sources. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --gene_phenotype
   - id: regulatory
     label: regulatory
     doc: |-
       Look for overlaps with regulatory regions. VEP can also report if a variant falls in a high information position within a transcription factor binding site. Output lines have a Feature type of RegulatoryFeature or MotifFeature. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --regulatory
   - id: cellType
     label: cellType
     doc: |-
       Report only regulatory regions that are found in the given cell type(s). Can be a single cell type or a comma-separated list. The functional type in each cell type is reported under CELL_TYPE in the output. To retrieve a list of cell types, use --cell_type list. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --cell_type
   - id: individual
     label: individual
     doc: |-
       Consider only alternate alleles present in the genotypes of the specified individual(s). May be a single individual, a comma-separated list or "all" to assess all individuals separately. Individual variant combinations homozygous for the given reference allele will not be reported. Each individual and variant combination is given on a separate line of output. Only works with VCF files containing individual genotype data; individual IDs are taken from column headers. Not used by default
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --individual
       itemSeparator: ','
   - id: phased
     label: phased
     doc: |-
       Force VCF genotypes to be interpreted as phased. For use with plugins that depend on phased data. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --phased
   - id: alleleNumber
     label: alleleNumber
     doc: |-
       Identify allele number from VCF input, where 1 = first ALT allele, 2 = second ALT allele etc. Useful when using --minimal Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --allele_number
   - id: showRefAllele
     label: showRefAllele
     doc: |-
       Adds the reference allele in the output. Mainly useful for the VEP "default" and tab-delimited output formats. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --show_ref_allele
   - id: totalLength
     label: totalLength
     doc: Give cDNA, CDS and protein positions as Position/Length. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --total_length
   - id: numbers
     label: numbers
     doc: |-
       Adds affected exon and intron numbering to to output. Format is Number/Total. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --numbers
   - id: noEscape
     label: noEscape
     doc: Don't URI escape HGVS strings. Default = escape
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_escape
   - id: keepCsq
     label: keepCsq
     doc: Don't overwrite existing CSQ entry in VCF INFO field. Overwrites by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --keep_csq
   - id: vcfInfoField
     label: vcfInfoField
     doc: |-
       Change the name of the INFO key that VEP write the consequences to in its VCF output. Use "ANN" for compatibility with other tools such as snpEff. Default: CSQ
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --vcf_info_field
   - id: terms
     label: terms
     doc: |-
       (-t) The type of consequence terms to output. The Ensembl terms are described here. The Sequence Ontology is a joint effort by genome annotation centres to standardise descriptions of biological sequences. Default = "SO"
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --terms
   - id: noHeaders
     label: noHeaders
     doc: Don't write header lines in output files. Default = add headers
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_headers
   - id: hgvs
     label: hgvs
     doc: |-
       Add HGVS nomenclature based on Ensembl stable identifiers to the output. Both coding and protein sequence names are added where appropriate. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. HGVS notations given on Ensembl identifiers are versioned. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --hgvs
   - id: hgvsg
     label: hgvsg
     doc: |-
       Add genomic HGVS nomenclature based on the input chromosome name. To generate HGVS identifiers when using --cache or --offline you must use a FASTA file and --fasta. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --hgvsg
   - id: shiftHgvs
     label: shiftHgvs
     doc: |-
       Enable or disable 3' shifting of HGVS notations. When enabled, this causes ambiguous insertions or deletions (typically in repetetive sequence tracts) to be "shifted" to their most 3' possible coordinates (relative to the transcript sequence and strand) before the HGVS notations are calculated; the flag HGVS_OFFSET is set to the number of bases by which the variant has shifted, relative to the input genomic coordinates. Disabling retains the original input coordinates of the variant. Default: 1 (shift)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --shift_hgvs
   - id: transcriptVersion
     label: transcriptVersion
     doc: Add version numbers to Ensembl transcript identifiers
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --transcript_version
   - id: protein
     label: protein
     doc: |-
       Add the Ensembl protein identifier to the output where appropriate. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --protein
   - id: symbol
     label: symbol
     doc: |-
       Adds the gene symbol (e.g. HGNC) (where available) to the output. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --symbol
   - id: ccds
     label: ccds
     doc: |-
       Adds the CCDS transcript identifer (where available) to the output. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --ccds
   - id: uniprot
     label: uniprot
     doc: |-
       Adds best match accessions for translated protein products from three UniProt-related databases (SWISSPROT, TREMBL and UniParc) to the output. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --uniprot
   - id: tsl
     label: tsl
     doc: |-
       Adds the transcript support level for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --tsl
   - id: appris
     label: appris
     doc: |-
       Adds the APPRIS isoform annotation for this transcript to the output. Not used by default. Note: Only available for human on the GRCh38 assembly
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --appris
   - id: canonical
     label: canonical
     doc: |-
       Adds a flag indicating if the transcript is the canonical transcript for the gene. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --canonical
   - id: mane
     label: mane
     doc: |-
       Adds a flag indicating if the transcript is the MANE Select transcript for the gene. Not used by default. Note: Only available for human on the GRCh38 assembly
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --mane
   - id: biotype
     label: biotype
     doc: Adds the biotype of the transcript or regulatory feature. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --biotype
   - id: domains
     label: domains
     doc: Adds names of overlapping protein domains to output. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --domains
   - id: xrefRefseq
     label: xrefRefseq
     doc: |-
       Output aligned RefSeq mRNA identifier for transcript. Not used by default. Note: The RefSeq and Ensembl transcripts aligned in this way MAY NOT, AND FREQUENTLY WILL NOT, match exactly in sequence, exon structure and protein product
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --xref_refseq
   - id: synonyms
     label: synonyms
     doc: |-
       Load a file of chromosome synonyms. File should be tab-delimited with the primary identifier in column 1 and the synonym in column 2. Synonyms allow different chromosome identifiers to be used in the input file and any annotation source (cache, database, GFF, custom file, FASTA file). Not used by default
     type:
     - File
     - 'null'
     inputBinding:
       prefix: --synonyms
   - id: checkExisting
     label: checkExisting
     doc: |-
       Checks for the existence of known variants that are co-located with your input. By default the alleles are compared and variants on an allele-specific basis - to compare only coordinates, use --no_check_alleles.

                   Some databases may contain variants with unknown (null) alleles and these are included by default; to exclude them use --exclude_null_alleles.

                   See this page for more details.

                   Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --check_existing
   - id: checkSvs
     label: checkSvs
     doc: |-
       Checks for the existence of structural variants that overlap your input. Currently requires database access. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --check_svs
   - id: clinSigAllele
     label: clinSigAllele
     doc: |-
       Return allele specific clinical significance. Setting this option to 0 will provide all known clinical significance values at the given locus. Default: 1 (Provide allele-specific annotations)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --clin_sig_allele
   - id: excludeNullAlleles
     label: excludeNullAlleles
     doc: |-
       Do not include variants with unknown alleles when checking for co-located variants. Our human database contains variants from HGMD and COSMIC for which the alleles are not publically available; by default these are included when using --check_existing, use this flag to exclude them. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude_null_alleles
   - id: noCheckAlleles
     label: noCheckAlleles
     doc: |-
       When checking for existing variants, by default VEP only reports a co-located variant if none of the input alleles are novel. For example, if your input variant has alleles A/G, and an existing co-located variant has alleles A/C, the co-located variant will not be reported.

                   Strand is also taken into account - in the same example, if the input variant has alleles T/G but on the negative strand, then the co-located variant will be reported since its alleles match the reverse complement of input variant.

                   Use this flag to disable this behaviour and compare using coordinates alone. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_check_alleles
   - id: af
     label: af
     doc: |-
       Add the global allele frequency (AF) from 1000 Genomes Phase 3 data for any known co-located variant to the output. For this and all --af_* flags, the frequency reported is for the input allele only, not necessarily the non-reference or derived allele. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --af
   - id: maxAf
     label: maxAf
     doc: |-
       Report the highest allele frequency observed in any population from 1000 genomes, ESP or gnomAD. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --max_af
   - id: af1kg
     label: af1kg
     doc: |-
       Add allele frequency from continental populations (AFR,AMR,EAS,EUR,SAS) of 1000 Genomes Phase 3 to the output. Must be used with --cache. Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --af_1kg
   - id: afEsp
     label: afEsp
     doc: |-
       Include allele frequency from NHLBI-ESP populations. Must be used with --cache. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --af_esp
   - id: afGnomad
     label: afGnomad
     doc: |-
       Include allele frequency from Genome Aggregation Database (gnomAD) exome populations. Note only data from the gnomAD exomes are included; to retrieve data from the additional genomes data set, see this guide. Must be used with --cache Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --af_gnomad
   - id: afExac
     label: afExac
     doc: |-
       Include allele frequency from ExAC project populations. Must be used with --cache. Not used by default. Note: ExAC data has been superceded by gnomAD. This flag remains for those wishing to use older cache versions containing ExAC data.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --af_exac
   - id: pubmed
     label: pubmed
     doc: |-
       Report Pubmed IDs for publications that cite existing variant. Must be used with --cache. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pubmed
   - id: failed
     label: failed
     doc: |-
       When checking for co-located variants, by default VEP will exclude variants that have been flagged as failed. Set this flag to include such variants. Default: 0 (exclude)
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --failed
   - id: gencodeBasic
     label: gencodeBasic
     doc: |-
       Limit your analysis to transcripts belonging to the GENCODE basic set. This set has fragmented or problematic transcripts removed. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --gencode_basic
   - id: excludePredicted
     label: excludePredicted
     doc: |-
       When using the RefSeq or merged cache, exclude predicted transcripts (i.e. those with identifiers beginning with "XM_" or "XR_").
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --exclude_predicted
   - id: transcriptFilter
     label: transcriptFilter
     doc: |-
       ADVANCED Filter transcripts according to any arbitrary set of rules. Uses similar notation to filter_vep.

                   You may filter on any key defined in the root of the transcript object; most commonly this will be ""stable_id"":

                   --transcript_filter ""stable_id match N[MR]_""
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --transcript_filter
   - id: checkRef
     label: checkRef
     doc: |-
       Force VEP to check the supplied reference allele against the sequence stored in the Ensembl Core database or supplied FASTA file. Lines that do not match are skipped. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --check_ref
   - id: lookupRef
     label: lookupRef
     doc: |-
       Force overwrite the supplied reference allele with the sequence stored in the Ensembl Core database or supplied FASTA file. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --lookup_ref
   - id: dontSkip
     label: dontSkip
     doc: |-
       Don't skip input variants that fail validation, e.g. those that fall on unrecognised sequences. Combining --check_ref with --dont_skip will add a CHECK_REF output field when the given reference does not match the underlying reference sequence.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --dont_skip
   - id: allowNonVariant
     label: allowNonVariant
     doc: |-
       When using VCF format as input and output, by default VEP will skip non-variant lines of input (where the ALT allele is null). Enabling this option the lines will be printed in the VCF output with no consequence data added.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --allow_non_variant
   - id: chr
     label: chr
     doc: |-
       Select a subset of chromosomes to analyse from your file. Any data not on this chromosome in the input will be skipped. The list can be comma separated, with "-" characters representing an interval. For example, to include chromosomes 1, 2, 3, 10 and X you could use --chr 1-3,10,X Not used by default
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --chr
       itemSeparator: ','
   - id: codingOnly
     label: codingOnly
     doc: |-
       Only return consequences that fall in the coding regions of transcripts. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --coding_only
   - id: noIntergenic
     label: noIntergenic
     doc: Do not include intergenic consequences in the output. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --no_intergenic
   - id: pick
     label: pick
     doc: |-
       Pick once line or block of consequence data per variant, including transcript-specific columns. Consequences are chosen according to the criteria described here, and the order the criteria are applied may be customised with --pick_order. This is the best method to use if you are interested only in one consequence per variant. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pick
   - id: pickAllele
     label: pickAllele
     doc: |-
       Like --pick, but chooses one line or block of consequence data per variant allele. Will only differ in behaviour from --pick when the input variant has multiple alternate alleles. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pick_allele
   - id: perGene
     label: perGene
     doc: |-
       Output only the most severe consequence per gene. The transcript selected is arbitrary if more than one has the same predicted consequence. Uses the same ranking system as --pick. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --per_gene
   - id: pickAlleleGene
     label: pickAlleleGene
     doc: |-
       Like --pick_allele, but chooses one line or block of consequence data per variant allele and gene combination. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --pick_allele_gene
   - id: flagPick
     label: flagPick
     doc: |-
       As per --pick, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --flag_pick
   - id: flagPickAllele
     label: flagPickAllele
     doc: |-
       As per --pick_allele, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --flag_pick_allele
   - id: flagPickAlleleGene
     label: flagPickAlleleGene
     doc: |-
       As per --pick_allele_gene, but adds the PICK flag to the chosen block of consequence data and retains others. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --flag_pick_allele_gene
   - id: pickOrder
     label: pickOrder
     doc: |-
       Customise the order of criteria (and the list of criteria) applied when choosing a block of annotation data with one of the following options: --pick, --pick_allele, --per_gene, --pick_allele_gene, --flag_pick, --flag_pick_allele, --flag_pick_allele_gene. See this page for the default order.
                   Valid criteria are: [ canonical appris tsl biotype ccds rank length mane ]. e.g.:

                   --pick --pick_order tsl,appris,rank
     type:
     - type: array
       items: string
     - 'null'
     inputBinding:
       prefix: --pick_order
       itemSeparator: ','
   - id: mostSevere
     label: mostSevere
     doc: |-
       Output only the most severe consequence per variant. Transcript-specific columns will be left blank. Consequence ranks are given in this table. To include regulatory consequences, use the --regulatory option in combination with this flag. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --most_severe
   - id: summary
     label: summary
     doc: |-
       Output only a comma-separated list of all observed consequences per variant. Transcript-specific columns will be left blank. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --summary
   - id: filterCommon
     label: filterCommon
     doc: |-
       Shortcut flag for the filters below - this will exclude variants that have a co-located existing variant with global AF > 0.01 (1%). May be modified using any of the following freq_* filters. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --filter_common
   - id: checkFrequency
     label: checkFrequency
     doc: |-
       Turns on frequency filtering. Use this to include or exclude variants based on the frequency of co-located existing variants in the Ensembl Variation database. You must also specify all of the --freq_* flags below. Frequencies used in filtering are added to the output under the FREQS key in the Extra field. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --check_frequency
   - id: freqPop
     label: freqPop
     doc: |-
       Name of the population to use in frequency filter. This must be one of the following: (1KG_ALL, 1KG_AFR, 1KG_AMR, 1KG_EAS, 1KG_EUR, 1KG_SAS, AA, EA, gnomAD, gnomAD_AFR, gnomAD_AMR, gnomAD_ASJ, gnomAD_EAS, gnomAD_FIN, gnomAD_NFE, gnomAD_OTH, gnomAD_SAS)
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --freq_pop
   - id: freqFreq
     label: freqFreq
     doc: Allele frequency to use for filtering. Must be a float value between 0 and
       1
     type:
     - float
     - 'null'
     inputBinding:
       prefix: --freq_freq
   - id: freqGtLt
     label: freqGtLt
     doc: |-
       Specify whether the frequency of the co-located variant must be greater than (gt) or less than (lt) the value specified with --freq_freq
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --freq_gt_lt
   - id: freqFilter
     label: freqFilter
     doc: |-
       Specify whether to exclude or include only variants that pass the frequency filter
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --freq_filter
   - id: caddReference
     label: caddReference
     type:
     - type: array
       items: File
     - 'null'
   - id: condelConfig
     label: condelConfig
     doc: "Directory containing CondelPlugin config, in format: '<dir>/condel_SP.conf'"
     type:
     - Directory
     - 'null'
   - id: dbnspReference
     label: dbnspReference
     doc: ''
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
   - id: dbsnpColumns
     label: dbsnpColumns
     type:
     - type: array
       items: string
     - 'null'
   - id: revelReference
     label: revelReference
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
   - id: custom1Reference
     label: custom1Reference
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
   - id: custom1Columns
     label: custom1Columns
     type:
     - type: array
       items: string
     - 'null'
   - id: custom2Reference
     label: custom2Reference
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
   - id: custom2Columns
     label: custom2Columns
     type:
     - type: array
       items: string
     - 'null'
   - id: database
     label: database
     doc: Enable VEP to use local or remote databases.
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --database
   - id: host
     label: host
     doc: |-
       Manually define the database host to connect to. Users in the US may find connection and transfer speeds quicker using our East coast mirror, useastdb.ensembl.org. Default = "ensembldb.ensembl.org"
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --host
   - id: user
     label: user
     doc: (-u) Manually define the database username. Default = "anonymous"
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --user
   - id: password
     label: password
     doc: (--pass) Manually define the database password. Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --password
   - id: port
     label: port
     doc: Manually define the database port. Default = 5306
     type:
     - int
     - 'null'
     inputBinding:
       prefix: --port
   - id: genomes
     label: genomes
     doc: |-
       Override the default connection settings with those for the Ensembl Genomes public MySQL server. Required when using any of the Ensembl Genomes species. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --genomes
   - id: isMultispecies
     label: isMultispecies
     doc: |-
       Some of the Ensembl Genomes databases (mainly bacteria and protists) are composed of a collection of close species. It updates the database connection settings (i.e. the database name) if the value is set to 1. Default: 0
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --is_multispecies
   - id: lrg
     label: lrg
     doc: |-
       Map input variants to LRG coordinates (or to chromosome coordinates if given in LRG coordinates), and provide consequences on both LRG and chromosomal transcripts. Not used by default
     type:
     - boolean
     - 'null'
     inputBinding:
       prefix: --lrg
   - id: dbVersion
     label: dbVersion
     doc: |-
       Force VEP to connect to a specific version of the Ensembl databases. Not recommended as there may be conflicts between software and database versions. Not used by default
     type:
     - string
     - 'null'
     inputBinding:
       prefix: --db_version
   - id: registry
     label: registry
     doc: |-
       Defining a registry file overwrites other connection settings and uses those found in the specified registry file to connect. Not used by default
     type:
     - string
     - 'null'
     default: generated
     inputBinding:
       prefix: --registry

   outputs:
   - id: out
     label: out
     type:
     - File
     - 'null'
     secondaryFiles:
     - pattern: .tbi
     outputBinding:
       glob: $(inputs.inputFile.basename.replace(/.vcf.gz$/, "")).vcf
       loadContents: false
   - id: out_stdout
     label: out_stdout
     type: stdout
   - id: out_stats
     label: out_stats
     type:
     - File
     - 'null'
     outputBinding:
       glob: |-
         $((inputs.statsFile ? inputs.statsFile : "generated" != null) ? inputs.statsFile ? inputs.statsFile : "generated" : "variant_effect_output.txt_summary.html")
       loadContents: false
   - id: out_warnings
     label: out_warnings
     type:
     - File
     - 'null'
     outputBinding:
       glob: generated-warning.txt
       loadContents: false
   stdout: _stdout
   stderr: _stderr

   baseCommand: vep
   arguments:
   - position: 0
     valueFrom: |-
       $((inputs.caddReference != null) ? ("--plugin CADD," + inputs.caddReference.join(",")) : "")
     shellQuote: false
   - position: 0
     valueFrom: |-
       $((inputs.condelConfig != null) ? "--plugin Condel,{condelconfig},b".replace(/\{condelconfig\}/g, inputs.condelConfig) : "")
     shellQuote: false
   - position: 0
     valueFrom: |-
       $(((inputs.dbnspReference != null) && (inputs.dbsnpColumns != null)) ? "--plugin dbNSFP,{ref},{cols}".replace(/\{ref\}/g, inputs.dbnspReference).replace(/\{cols\}/g, inputs.dbsnpColumns.join(",")) : "")
     shellQuote: false
   - position: 0
     valueFrom: |-
       $((inputs.revelReference != null) ? "--plugin REVEL,{ref}".replace(/\{ref\}/g, inputs.revelReference) : "")
     shellQuote: false
   - position: 0
     valueFrom: |-
       $(((inputs.custom1Reference != null) && (inputs.custom1Columns != null)) ? "--custom {ref},{cols}".replace(/\{ref\}/g, inputs.custom1Reference).replace(/\{cols\}/g, inputs.custom1Columns.join(",")) : "")
     shellQuote: false
   - position: 0
     valueFrom: |-
       $(((inputs.custom2Reference != null) && (inputs.custom2Columns != null)) ? "--custom {ref},{cols}".replace(/\{ref\}/g, inputs.custom2Reference).replace(/\{cols\}/g, inputs.custom2Columns.join(",")) : "")
     shellQuote: false

   hints:
   - class: ToolTimeLimit
     timelimit: |-
       $([inputs.runtime_seconds, 86400].filter(function (inner) { return inner != null })[0])
   id: vep_database


