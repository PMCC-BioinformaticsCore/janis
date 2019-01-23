
strelka
=======
*bioinformatics*

Documentation
-------------

URL
******
`https://github.com/Illumina/strelka <https://github.com/Illumina/strelka>`_

Docstring
*********
Strelka2 is a fast and accurate small variant caller optimized for analysis of germline variation 
    in small cohorts and somatic variation in tumor/normal sample pairs. The germline caller employs 
    an efficient tiered haplotype model to improve accuracy and provide read-backed phasing, adaptively 
    selecting between assembly and a faster alignment-based haplotyping approach at each variant locus. 
    The germline caller also analyzes input sequencing data using a mixture-model indel error estimation 
    method to improve robustness to indel noise. The somatic calling model improves on the original 
    Strelka method for liquid and late-stage tumor analysis by accounting for possible tumor cell 
    contamination in the normal sample. A final empirical variant re-scoring step using random forest 
    models trained on various call quality features has been added to both callers to further improve precision.
    
    Compared with submissions to the recent PrecisonFDA Consistency and Truth challenges, the average 
    indel F-score for Strelka2 running in its default configuration is 3.1% and 0.08% higher, respectively, 
    than the best challenge submissions. Runtime on a 28-core server is ~40 minutes for 40x WGS germline 
    analysis and ~3 hours for a 110x/40x WGS tumor-normal somatic analysis
    
    Strelka accepts input read mappings from BAM or CRAM files, and optionally candidate and/or forced-call 
    alleles from VCF. It reports all small variant predictions in VCF 4.1 format. Germline variant 
    reporting uses the gVCF conventions to represent both variant and reference call confidence. 
    For best somatic indel performance, Strelka is designed to be run with the Manta structural variant 
    and indel caller, which provides additional indel candidates up to a given maxiumum indel size 
    (49 by default). By design, Manta and Strelka run together with default settings provide complete 
    coverage over all indel sizes (in additional to SVs and SNVs). 
    
    See the user guide for a full description of capabilities and limitations

Outputs
-------
============  ==========  ===========================================================================================================================================================================================================================================
name          type        documentation
============  ==========  ===========================================================================================================================================================================================================================================
directory     Directory
configPickle  File
script        File
stats         tsv         A tab-delimited report of various internal statistics from the variant calling process: Runtime information accumulated for each genome segment, excluding auxiliary steps such as BAM indexing and vcf merging. Indel candidacy statistics
variants      vcf-gz-tbi  Primary variant inferences are provided as a series of VCF 4.1 files
genome        vcf-gz-tbi
============  ==========  ===========================================================================================================================================================================================================================================

Inputs
------
========================  =================  ================  ==========  =========================================================================================================================================================================================================
name                      type               prefix              position  documentation
========================  =================  ================  ==========  =========================================================================================================================================================================================================
bam                       BamPair            --bam                      1  Sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [required] (no default)
reference                 FastaWithDict      --referenceFasta           1  samtools-indexed reference fasta file [required]
relativeStrelkaDirectory  Optional<String>   --runDir                   1  Name of directory to be created where all workflow scripts and output will be written. Each analysis requires a separate directory.
version                   Optional<Boolean>  --version                  3  show program's version number and exit
help                      Optional<Boolean>  --help                     3  (-h) show this help message and exit
mode                      Optional<String>   --mode                     3  (-m MODE)  select run mode (local|sge)
queue                     Optional<String>   --queue                    3  (-q QUEUE) specify scheduler queue name
jobs                      Optional<String>   --jobs                     3  (-j JOBS)  number of jobs, must be an integer or 'unlimited' (default: Estimate total cores on this node for local mode, 128 for sge mode)
memGb                     Optional<String>   --memGb                    3  (-g MEMGB) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local mode, 'unlimited' for sge mode)
dryRun                    Optional<Boolean>  --dryRun                   3  dryRun (-d,) workflow code without actually running command-tasks
quiet                     Optional<Boolean>  --quiet                    3  Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
mailTo                    Optional<String>   --mailTo                   3  (-e) send email notification of job completion status to this address (may be provided multiple times for more than one email address)
========================  =================  ================  ==========  =========================================================================================================================================================================================================


*This page was automatically generated*
