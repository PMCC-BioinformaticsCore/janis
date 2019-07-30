:orphan:


Manta
=============

Description
-------------

Tool identifier: ``manta``

Tool path: ``janis_bioinformatics.tools.illumina.manta.manta import Manta_1_4_0``

Version: 1.4.0

Container: ``michaelfranklin/manta:1.4.0``

Versions
*********

- `1.5.0 <manta_1.5.0.html>`_
- 1.4.0 (current)

Documentation
-------------

URL
******
`https://github.com/Illumina/manta <https://github.com/Illumina/manta>`_

Tool documentation
******************
Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads. 
It is optimized for analysis of germline variation in small sets of individuals and somatic 
variation in tumor/normal sample pairs. Manta discovers, assembles and scores large-scale SVs, 
medium-sized indels and large insertions within a single efficient workflow. The method is 
designed for rapid analysis on standard compute hardware: NA12878 at 50x genomic coverage is 
analyzed in less than 20 minutes on a 20 core server, and most WGS tumor/normal analyses 
can be completed within 2 hours. Manta combines paired and split-read evidence during SV 
discovery and scoring to improve accuracy, but does not require split-reads or successful 
breakpoint assemblies to report a variant in cases where there is strong evidence otherwise. 

It provides scoring models for germline variants in small sets of diploid samples and somatic 
variants in matched tumor/normal sample pairs. There is experimental support for analysis of 
unmatched tumor samples as well. Manta accepts input read mappings from BAM or CRAM files and 
reports all SV and indel inferences in VCF 4.1 format. See the user guide for a full description 
of capabilities and limitations.

Outputs
-------
==========================  ==========  ===============
name                        type        documentation
==========================  ==========  ===============
python                      File
pickle                      File
candidateSV                 vcf-gz-tbi
candidateSmallIndels        vcf-gz-tbi
diploidSV                   vcf-gz-tbi
alignmentStatsSummary       File
svCandidateGenerationStats  tsv
svLocusGraphStats           tsv
==========================  ==========  ===============

Inputs
------
Find the inputs below

Required inputs
***************

=========  =============  ================  ==========  ===============================================================================================================================================================================
name       type           prefix              position  documentation
=========  =============  ================  ==========  ===============================================================================================================================================================================
bam        BamPair        --bam                      1  FILE Normal sample BAM or CRAM file. May be specified more than once, multiple inputs will be treated as each BAM file representing a different sample. [optional] (no default)
reference  FastaWithDict  --referenceFasta           1  samtools-indexed reference fasta file [required]
=========  =============  ================  ==========  ===============================================================================================================================================================================

Optional inputs
***************

==============  ==================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================
name            type                prefix              position  documentation
==============  ==================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================
config          Optional<File>      --config                   1  provide a configuration file to override defaults in global config file (/opt/conda/share/manta-1.2.1-0/bin/configManta.py.ini)
runDir          Optional<Filename>  --runDir                   1  Run script and run output will be written to this directory [required] (default: MantaWorkflow)
tumorBam        Optional<BamPair>   --tumorBam                 1  Tumor sample BAM or CRAM file. Only up to one tumor bam file accepted. [optional=null]
exome           Optional<Boolean>   --exome                    1  Set options for WES input: turn off depth filters
rna             Optional<BAM>       --rna                      1  Set options for RNA-Seq input. Must specify exactly one bam input file
unstrandedRNA   Optional<File>      --unstrandedRNA            1  Set if RNA-Seq input is unstranded: Allows splice-junctions on either strand
outputContig    Optional<File>      --outputContig             1  Output assembled contig sequences in VCF file
callRegions     Optional<BedTABIX>  --callRegions              1  Optionally provide a bgzip-compressed/tabix-indexed BED file containing the set of regions to call. No VCF output will be provided outside of these regions. The full genome will still be used to estimate statistics from the input (such as expected depth per chromosome). Only one BED file may be specified. (default: call the entire genome)
mode            Optional<String>    --mode                     3  (-m) select run mode (local|sge)
quiet           Optional<Boolean>   --quiet                    3  Don't write any log output to stderr (but still write to workspace/pyflow.data/logs/pyflow_log.txt)
queue           Optional<String>    --queue                    3  (-q) specify scheduler queue name
memgb           Optional<Integer>   --memGb                    3  (-g) gigabytes of memory available to run workflow -- only meaningful in local mode, must be an integer (default: Estimate the total memory for this node for local  mode, 'unlimited' for sge mode)
maxTaskRuntime  Optional<String>    --maxTaskRuntime           3  (format: hh:mm:ss) Specify scheduler max runtime per task, argument is provided to the 'h_rt' resource limit if using SGE (no default)
==============  ==================  ================  ==========  ====================================================================================================================================================================================================================================================================================================================================================


Metadata
********

Author: Michael Franklin


*Manta was last updated on 2019-02-19*.
*This page was automatically generated on 2019-07-30*.
