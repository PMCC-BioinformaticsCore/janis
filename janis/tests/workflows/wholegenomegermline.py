import unittest

from janis import Input, String, Step, Directory, Workflow, Array, Output
from bioinformatics import FastaWithDict
from bioinformatics import Fastq
from bioinformatics import VcfIdx, TabixIdx
from bioinformatics import BcfToolsNormLatest as BcfToolsNorm
from bioinformatics import AlignSortedBam
from bioinformatics import Gatk4ApplyBqsrLatest as Gatk4ApplyBqsr
from bioinformatics import \
    Gatk4BaseRecalibratorLatest as Gatk4BaseRecalibrator
from bioinformatics import \
    Gatk4HaplotypeCallerLatest as Gatk4HaplotypeCaller
from bioinformatics import Gatk4MarkDuplicatesLatest as Gatk4MarkDuplicates
from bioinformatics import Gatk4MergeSamFilesLatest as Gatk4MergeSamFiles
from bioinformatics import PerformanceValidator_1_2_1


class TestGermlinePipeline(unittest.TestCase):

    def test_workflow(self):

        w = Workflow("whole_genome_germline")

        reference = Input("reference", FastaWithDict())
        fastqInputs = Input("inputs", Array(Fastq()))

        s1_inp_header = Input("read_group_header_line", String())
        s4_inp_SNPS_dbSNP = Input("SNPS_dbSNP", VcfIdx())
        s4_inp_SNPS_1000GP = Input("SNPS_1000GP", TabixIdx())
        s4_inp_OMNI = Input("OMNI", TabixIdx())
        s4_inp_HAPMAP = Input("HAPMAP", TabixIdx())
        validator_truth = Input("TRUTH_VCF", VcfIdx())
        validator_intervals = Input("INTERVALS", Array(VcfIdx()))

        inp_tmpdir = Input("tmpdir", Directory())

        s1_sw = Step("sw_bwa_st_sort", AlignSortedBam())
        s2_mergeSamFiles = Step("s2_mergeSamFiles", Gatk4MergeSamFiles())
        s3_markDuplicates = Step("s3_markDuplicates", Gatk4MarkDuplicates())
        s4_baseRecal = Step("s4_baseRecal", Gatk4BaseRecalibrator())
        s5_applyBqsr = Step("s5_applyBqsr", Gatk4ApplyBqsr())
        s6_haploy = Step("s6_haploy", Gatk4HaplotypeCaller())
        s7_bcfNorm = Step("s7_bcfNorm", BcfToolsNorm())
        s8_validator = Step("s8_validator", PerformanceValidator_1_2_1())

        # step1
        w.add_edge(fastqInputs, s1_sw.fastq)
        w.add_edges([
            (reference, s1_sw.reference),
            (s1_inp_header, s1_sw.read_group_header_line),
            (inp_tmpdir, s1_sw.tmpdir)
        ])

        #step2
        w.add_edge(s1_sw.o3_sortsam, s2_mergeSamFiles.input)
        w.add_edge(inp_tmpdir, s2_mergeSamFiles.tmpDir)
        w.add_default_value(s2_mergeSamFiles.useThreading, True)
        w.add_default_value(s2_mergeSamFiles.createIndex, True)
        w.add_default_value(s2_mergeSamFiles.maxRecordsInRam, 5000000)
        w.add_default_value(s2_mergeSamFiles.validationStringency, "SILENT")




        # step3
        w.add_edge(s2_mergeSamFiles, s3_markDuplicates)
        w.add_edge(inp_tmpdir, s3_markDuplicates.tmpDir)
        w.add_default_value(s3_markDuplicates.createIndex, True)
        w.add_default_value(s3_markDuplicates.maxRecordsInRam, 5000000)

        #step4 - baserecal
        w.add_edges([
            (s3_markDuplicates, s4_baseRecal),
            (reference, s4_baseRecal.reference),
            (s4_inp_SNPS_dbSNP, s4_baseRecal.knownSites),
            (s4_inp_SNPS_1000GP, s4_baseRecal.knownSites),
            (s4_inp_OMNI, s4_baseRecal.knownSites),
            (s4_inp_HAPMAP, s4_baseRecal.knownSites),
            (inp_tmpdir, s4_baseRecal.tmpDir)
        ])

        # step5 - apply bqsr
        w.add_edges([
            (s3_markDuplicates.output, s5_applyBqsr),
            (s4_baseRecal, s5_applyBqsr),
            (reference, s5_applyBqsr.reference),
            (inp_tmpdir, s5_applyBqsr.tmpDir)
        ])

        # step6 - haplotype caller
        w.add_edges([
            (s5_applyBqsr, s6_haploy),
            (reference, s6_haploy),
            (s4_inp_SNPS_dbSNP, s6_haploy)
        ])

        # step7 - BcfToolsNorm
        w.add_edges([
            (reference, s7_bcfNorm),
            (s6_haploy, s7_bcfNorm.input)
            # (s7_cheat_inp, s7_bcfNorm.input)
        ])

        # step8 - validator

        w.add_edges([
            (s7_bcfNorm, s8_validator),
            (validator_truth, s8_validator.truth),
            (validator_intervals, s8_validator.intervals)
        ])


        # Outputs

        w.add_edges([
            (s1_sw.o1_bwa, Output("sw_bwa")),
            (s1_sw.o2_samtools, Output("sw_samtools")),
            (s1_sw.o3_sortsam, Output("sw_sortsam")),
            (s2_mergeSamFiles, Output("o4_merged")),
            (s3_markDuplicates.output, Output("o5_marked_output")),
            (s3_markDuplicates.metrics, Output("o5_marked_metrics")),
            (s4_baseRecal, Output("o6_recal")),
            (s5_applyBqsr, Output("o7_bqsr")),
            (s6_haploy.output, Output("o8_halpo")),
            (s7_bcfNorm, Output("o9_bcfnorm")),
            (s8_validator.summaryMetrics, Output("o12_concord_summary")),
            (s8_validator.detailMetrics, Output("o12_concord_detail")),
            (s8_validator.contingencyMetrics, Output("o12_concord_contig"))
        ])

        # w.draw_graph()
        w.dump_cwl(to_disk=True, with_docker=True)
        w.dump_wdl(to_disk=True, with_docker=False)



original_bash = """
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --partition=debug
#SBATCH --mem=64G
#SBATCH --time=10-00:00:00
#SBATCH --mail-user=jiaan.yu@petermac.org
#SBATCH --mail-type=ALL
#SBATCH --output=pipeline-%j.out
#SBATCH --job-name="WGSpipeline"


echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

STARTTIME=$(date +%s)


module load bwa/0.7.17
module load samtools/1.9
module load gatk/4.0.10.0
module load bcftools/1.9
module load igvtools

## Input parameters: R1 of fastq files

READ_GROUP_HEADER_LINE='@RG\tID:NA12878\tSM:NA12878\tLB:NA12878\tPL:ILLUMINA'
REFERENCE="/bioinf_core/Proj/hg38_testing/Resources/Gatk_Resource_Bundle_hg38/hg38_contigs_renamed/Homo_sapiens_assembly38.fasta"
SNPS_dbSNP='/bioinf_core/Proj/hg38_testing/Resources/Gatk_Resource_Bundle_hg38/Homo_sapiens_assembly38.dbsnp138.vcf'
SNPS_1000GP='/bioinf_core/Proj/hg38_testing/Resources/Gatk_Resource_Bundle_hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz'
OMNI='/bioinf_core/Proj/hg38_testing/Resources/Gatk_Resource_Bundle_hg38/1000G_omni2.5.hg38.vcf.gz'
HAPMAP='/bioinf_core/Proj/hg38_testing/Resources/Gatk_Resource_Bundle_hg38/hapmap_3.3.hg38.vcf.gz'
TRUTH_VCF=${1/_R1.fastq/.vcf}
INTERVAL_LIST=${1/_R1.fastq/.interval_list}

## Fill in the output file here
SAM='test3/test.sam'


bwa mem -R $READ_GROUP_HEADER_LINE \
 -M \
 -t 36 \
 $REFERENCE \
 $1 ${1/R1/R2} >  $SAM


samtools view -S -h -b $SAM > ${SAM/sam/bam}


gatk SortSam \
 -I=$SAM \
 -O=${SAM/sam/sorted.bam} \
 --SORT_ORDER=coordinate \
 --VALIDATION_STRINGENCY=SILENT \
 --CREATE_INDEX=true \
 --MAX_RECORDS_IN_RAM=5000000 \
 --TMP_DIR='/researchers/jiaan.yu/WGS_pipeline/germline/tmp'


gatk MergeSamFiles \
 -I=${SAM/sam/sorted.bam} \
 -O=${SAM/sam/sorted.merged.bam} \
 --USE_THREADING=true \
 --CREATE_INDEX=true \
 --MAX_RECORDS_IN_RAM=5000000 \
 --VALIDATION_STRINGENCY=SILENT \
 --TMP_DIR='/researchers/jiaan.yu/WGS_pipeline/germline/tmp'


gatk MarkDuplicates \
 -I=${SAM/sam/sorted.merged.bam} \
 -O=${SAM/sam/sorted.merged.markdups.bam} \
 --CREATE_INDEX=true \
 --METRICS_FILE=${SAM/sam/metrics.txt} \
 --MAX_RECORDS_IN_RAM=5000000 \
 --TMP_DIR='/researchers/jiaan.yu/WGS_pipeline/germline/tmp' 


gatk BaseRecalibrator \
 -I=${SAM/sam/sorted.merged.markdups.bam} \
 -O=${SAM/sam/recal.table} \
 -R=$REFERENCE \
 --known-sites=$SNPS_dbSNP \
 --known-sites=$SNPS_1000GP \
 --known-sites=$OMNI \
 --known-sites=$HAPMAP \
 --tmp-dir='/researchers/jiaan.yu/WGS_pipeline/germline/tmp'


gatk ApplyBQSR \
 -R=$REFERENCE \
 -I=${SAM/sam/sorted.merged.markdups.bam} \
 --bqsr-recal-file=${SAM/sam/recal.table} \
 -O=${SAM/sam/sorted.merged.markdups.recal.bam} \
 --tmp-dir='/researchers/jiaan.yu/WGS_pipeline/germline/tmp'


gatk HaplotypeCaller \
 -I=${SAM/sam/sorted.merged.markdups.recal.bam} \
 -R=$REFERENCE \
 -O=${SAM/sam/hap.vcf}
 -D=$SNPS_dbSNP


# The normalised vcf for future steps
bcftools norm -m -both -f $REFERENCE ${SAM/sam/hap.vcf} -o ${SAM/sam/hap.norm.vcf}

# The compressed and indexed vcf for the sake for the validation
bgzip -c ${SAM/sam/hap.norm.vcf} > ${SAM/sam/hap.norm.vcf.gz}
tabix -p vcf ${SAM/sam/hap.norm.vcf.gz}
gatk GenotypeConcordance \
 --TRUTH_VCF $TRUTH_VCF \
 --CALL_VCF ${SAM/sam/hap.norm.vcf.gz} \
 --OUTPUT ${SAM/sam/hap} \
 --MISSING_SITES_HOM_REF true \
 --INTERVALS $INTERVAL_LIST


ENDTIME=$(date +%s)

echo "Time to complete: $(($ENDTIME - $STARTTIME))"
"""
