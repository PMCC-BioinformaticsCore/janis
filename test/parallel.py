import examples.unix_commands
import examples.bio_informatics
from pipeline_definition.utils.logger import Logger, LogLevel


bwa_ref = ".bai, amb, ann, bwt, pac, sa, fai"
ref = ".fasta, .amb ann bwt pac sa fai ^dict"
dbsnp = ".vcf .vcf.ifx"

"""
Tools to implement:
    bwa-mem
    gather
    gatk-base-recalibrator
    gatk-haplotype
    gatk-mutect2
    gatk-printreads
    picard-markdups
    picard-sortsam
    samtools-view
"""


_yml = """
inputs:
    bwa_threads: 8
    bwa_tumor_reads:
        - type: File
          path: /Users/franklinmichael/source/seqliner-cwl/data/Vespa10-2029T-REDUCED_ACTTGA_L007_R1_001.fastq.gz 
        - type: File
          path: /Users/franklinmichael/source/seqliner-cwl/data/Vespa10-2029T-REDUCED_ACTTGA_L007_R2_001.fastq.gz 
    

    bwa_ref:
        type: RefFasta
        base: path/to/bwa/ref.fasta
        
    sam_tumor_output_name: string
    sam_normal_output_name: string
        
steps:

    # ====== TUMOUR ====== #

    bwa-mem-tumor:
        tool: bwa-mem
        vcfIdx: bwa_tumor_output_filename
        readGroup: bwa_readGroup
        reads: bwa_tumor_reads 
        reference: bwa_reference
        threads: bwa_threads
        
    sam-view-tumor:
        tool: samtools-view
        input: bwa-mem-tumor/out
        output_name: sam_tumor_output_name
        
    picard-sortsam-tumor:
        tool: picard-sortsam
        inputFileName_sortSam: sam-view-tumor/output
        outputFileName_sortSam: sortsam_tumor_outputFileName_sortSam
        tmpdir: sortsam_tmpdir
        validation_stringency: sortsam_validation_stringency
        maxRecordsInRam: sortsam_maxRecordsInRam
        
    picard-markdup-tumor:
        tool: picard-markdups
        inputFileName_markDups:
            source: picard-sortsam-tumor/output
            valueFrom: ${return [ self ];}
            
        outputFileName_markDups: markdups_tumor_outputFileName_markDups
        metricsFile: markdups_tumor_metricsFile
        picard_markdup_tmpdir: markdups_picard_markdup_tmpdir
        maxRecordsInRam: markdups_maxRecordsInRam
        validation_stringency: markdups_validation_stringency
        
    gather-md-tumor-index:
        tool: gather
        bamFile:
            source: picard-markdup-tumor/markDups_output
        bamIndex:
            source: picard-markdup-tumor/markDups_index_output 
            
    gatk-recalibrator-tumor:
        tool: gatk-base-recalibrator
        inputBam_BaseRecalibrator: gather-md-tumor-index/output
        outputfile_BaseRecalibrator: recalibrator_tumor_outputfile_BaseRecalibrator
        reference: recalibrator_reference
        known: recalibrator_known
        bedFile: recalibrator_bedFile
        
    gatk-printreads-tumor:
        tool: gatk-printreads
        inputBam_printReads: gather-md-tumor-index/output
        input_baseRecalibrator: gatk-recalibrator-tumor/output
        reference: printreads_reference
        outputfile_printReads: printreads_tumor_outputfile_printReads
        bedFile: printreads_bedFile
        # out: [o_printReads, printreads_idx_output]
        
    gather-pr-tumor-index:
        tool: gather
        bamFile:
            source: gatk-printreads-tumor/output_printReads
        bamIndex:
            source: gatk-printreads-tumor/printreads_index_output

    gatk-haplotypecaller-tumor:
        tool: gatk-haplotype
        inputBam_HaplotypeCaller: gather-pr-tumor-index/output
        reference: haplotype_reference
        outputfile_HaplotypeCaller: haplotype_tumor_outputfile_HaplotypeCaller
        dbsnp: haplotype_dbsnp
        threads: haplotype_threads
        emitRefConfidence: haplotype_emitRefConfidence
        bedFile: haplotype_bedFile
        bamOutput: haplotype_tumor_bamOutput
        
    # ====== NORMAL ====== #

    bwa-mem-normal:
        tool: bwa-mem
        output_filename: bwa_normal_output_filename
        readGroup: bwa_readGroup
        reads: bwa_normal_reads 
        reference: bwa_reference
        threads: bwa_threads

    sam-view-normal:
        tool: samtools-view
        input: bwa-mem-normal/output
        output_name: sam_normal_output_name
        
    picard-sortsam-normal:
        tool: picard-sortsam
        inputFileName_sortSam: sam-view-normal/output
        outputFileName_sortSam: sortsam_normal_outputFileName_sortSam
        tmpdir: sortsam_tmpdir
        validation_stringency: sortsam_validation_stringency
        maxRecordsInRam: sortsam_maxRecordsInRam
        
    picard-markdup-normal:
        tool: picard-markdups
        inputFileName_markDups: 
            source: picard-sortsam-normal/output
            valueFrom: ${return [ self ];}
        outputFileName_markDups: markdups_normal_outputFileName_markDups
        metricsFile: markdups_normal_metricsFile
        picard_markdup_tmpdir: markdups_picard_markdup_tmpdir
        maxRecordsInRam: markdups_maxRecordsInRam
        validation_stringency: markdups_validation_stringency

    gather-md-normal-index:
        tool: gather
        bamFile:
            source: picard-markdup-normal/markDups_output
        bamIndex:
            source: picard-markdup-normal/markDups_index_output
    
    gatk-recalibrator-normal:
        tool: gatk-base-recalibrator
        inputBam_BaseRecalibrator: gather-md-normal-index/output
        outputfile_BaseRecalibrator: recalibrator_normal_outputfile_BaseRecalibrator
        reference: recalibrator_reference
        known: recalibrator_known
        bedFile: recalibrator_bedFile
        
    gatk-printreads-normal:
        tool: gatk-printreads
        inputBam_printReads: gather-md-normal-index/output
        input_baseRecalibrator: gatk-recalibrator-normal/output
        reference: printreads_reference
        outputfile_printReads: printreads_normal_outputfile_printReads
        bedFile: printreads_bedFile
        # out: [output_printReads, printreads_index_output]
        
    gather-pr-normal-index:
        tool: gather
        bamFile: gatk-printreads-normal/output_printReads
        bamIndex: gatk-printreads-normal/printreads_index_output
        
    gatk-haplotypecaller-normal:
        tool: gatk-haplotype
        inputBam_HaplotypeCaller: gather-pr-normal-index/combined
        reference: haplotype_reference
        outputfile_HaplotypeCaller: haplotype_normal_outputfile_HaplotypeCaller
        dbsnp: haplotype_dbsnp
        threads: haplotype_threads
        emitRefConfidence: haplotype_emitRefConfidence
        bedFile: haplotype_bedFile
        bamOutput: haplotype_normal_bamOutput
    
    #=========================    

    gatk-mutect:
        tool: gatk-mutect2
        inputBam_tumor: gather-pr-tumor-index/output
        inputBam_normal: gather-pr-normal-index/output
        reference: mutect_reference
        bedFile: mutect_bedFile
        outputfile_name: mutect_outputfile_name
        dbsnp: mutect_dbsnp
        cosmic: mutect_cosmic
    
    #==================    
        
        
        
        
outputs:
    o1: bwa-mem-tumor/bwa_output
    o2: bwa-mem-normal/bwa_output
    o3: sam-view-tumor/sam_output
    o4: sam-view-normal/sam_output
    o5: picard-sortsam-tumor/sortSam_output
    o6: picard-sortsam-tumor/sortSam_output_indexes
    o7: picard-sortsam-tumor/sortSam_output
    o8: picard-sortsam-tumor/sortSam_output_indexes
    o9: picard-markdup-tumor/markDups_output
    o10: picard-markdup-tumor/markDups_metric_output
    o11: picard-markdup-tumor/markDups_index_output
    o12: picard-markdup-normal/markDups_output
    o13: picard-markdup-normal/markDups_metric_output
    o14: picard-markdup-normal/markDups_index_output
    o15: gatk-recalibrator-tumor/output_baseRecalibrator
    o16: gatk-recalibrator-normal/output_baseRecalibrator
    o17: gatk-printreads-tumor/output_printReads
    o18: gatk-printreads-tumor/printreads_index_output
    o19: gatk-printreads-normal/output_printReads
    o20: gatk-printreads-normal/printreads_index_output
    o21: gatk-haplotypecaller-tumor/output_HaplotypeCaller
    o22: gatk-haplotypecaller-tumor/bamOut_HaplotypeCaller
    o23: gatk-haplotypecaller-normal/output_HaplotypeCaller
    o24: gatk-haplotypecaller-normal/bamOut_HaplotypeCaller
    o25: gatk-mutect/output_mutect2
    o26: gather-md-normal-index/combined


"""

class ParallelPipeline():

    @staticmethod
    def test_simple():
        from pipeline_definition.pipeline_translator import PipelineTranslator
        Logger.set_write_location("/Users/franklinmichael/source/wehi-pipeline-definition/parallel.log")
        translator = PipelineTranslator()
        translation = translator.translate_string(_yml)
        Logger.close_file()
        print(translation)


if __name__ == '__main__':
    ParallelPipeline.test_simple()
