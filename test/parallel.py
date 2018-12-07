from utils.logger import Logger

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
hints:
    environment: person
    
inputs:
    bwa_threads: 8
    bwa_readGroup: string
    
    bwa_tumor_output_filename: string
    bwa_tumor_reads:
        - type: File
          path: /Users/franklinmichael/source/seqliner-cwl/data/Vespa10-2029T-REDUCED_ACTTGA_L007_R1_001.fastq.gz 
        - type: File
          path: /Users/franklinmichael/source/seqliner-cwl/data/Vespa10-2029T-REDUCED_ACTTGA_L007_R2_001.fastq.gz 
    
    bwa_normal_output_filename: string
    
    bwa_normal_reads:
        - type: File
          path: path/here
          
    bwa_ref:
        type: FastaRef
        base: path/to/bwa/ref.fasta
        
    sam_tumor_output_name: string
    sam_normal_output_name: string
    
    sortsam_tmpdir: string
    sortsam_validation_stringency: string
    sortsam_maxRecordsInRam: int
    sortsam_tumor_outputFileName_sortSam: string
    sortsam_normal_outputFileName_sortSam: string
    
    markdups_picard_markdup_tmpdir: string
    markdups_maxRecordsInRam: int
    markdups_validation_stringency: string
    markdups_tumor_outputFileName_markDups: string
    markdups_tumor_metricsFile: string
    markdups_normal_outputFileName_markDups: string
    markdups_normal_metricsFile: string
    
    recalibrator_reference:
        - type: FastaRef
          path: "path/here"
    
    recalibrator_known:
        - type: VcfIdx
          path: "path/here"
          
    recalibrator_bedFile: File
    recalibrator_tumor_outputfile_BaseRecalibrator: string
    recalibrator_normal_outputfile_BaseRecalibrator: string
    
    printreads_reference:
        - type: FastaRef
          path: "path/here"
          
    printreads_bedFile: File
    printreads_tumor_outputfile_printReads: string
    printreads_normal_outputfile_printReads: string
    
    haplotype_reference:
        - type: FastaRef
          path: "path/here"
    haplotype_dbsnp: 
        type: VcfIdx
          
    haplotype_threads: 5
    haplotype_emitRefConfidence: string
    haplotype_bedFile: 
        type: File
        path: something/here
        
    haplotype_tumor_bamOutput: string
    haplotype_normal_bamOutput: string

    haplotype_tumor_outputfile_HaplotypeCaller: string
    haplotype_normal_fOutput: string
    haplotype_normal_outputfile_HaplotypeCaller: string
    mutect_reference:
        type: FastaRef
        path: here
        
    mutect_bedFile: File
    mutect_outputfile_name: string
    mutect_dbsnp:
        - type: VcfIdx
          path: path/to/file
    mutect_cosmic:
        - type: VcfIdx
          path: path/to/file
    
steps:

    # ====== TUMOUR ====== #

    bwa-mem-tumor:
        tool: bwa-mem
        outputFilename: bwa_tumor_output_filename
        readGroup: bwa_readGroup
        reads: bwa_tumor_reads 
        reference: bwa_ref
        threads: bwa_threads
        
    sam-view-tumor:
        tool: samtools-view
        input: bwa-mem-tumor/out
        outputName: sam_tumor_output_name
        
    picard-sortsam-tumor:
        tool: picard-sortsam
        inputFileName_sortSam: sam-view-tumor/out
        outputFileName_sortSam: sortsam_tumor_outputFileName_sortSam
        tmpdir: sortsam_tmpdir
        validation_stringency: sortsam_validation_stringency
        maxRecordsInRam: sortsam_maxRecordsInRam
        
    picard-markdup-tumor:
        tool: picard-markdups
        inputFileName_markDups: picard-sortsam-tumor/out
           # source: picard-sortsam-tumor/out
           # valueFrom: ${return [ self ];}
            
        outputFileName_markDups: markdups_tumor_outputFileName_markDups
        metricsFile: markdups_tumor_metricsFile
        picard_markdup_tmpdir: markdups_picard_markdup_tmpdir
        maxRecordsInRam: markdups_maxRecordsInRam
        validation_stringency: markdups_validation_stringency
        hints: 
            
        
    gather-md-tumor-index:
        tool: gather
        bamFile: picard-markdup-tumor/out
        bamIndex: picard-markdup-tumor/out_idx 
            
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
        input_baseRecalibrator: gatk-recalibrator-tumor/out
        reference: printreads_reference
        outputfile_printReads: printreads_tumor_outputfile_printReads
        bedFile: printreads_bedFile
        # out: [o_printReads, printreads_idx_output]
        
    gather-pr-tumor-index:
        tool: gather
        bamFile: gatk-printreads-tumor/out
        bamIndex: gatk-printreads-tumor/out_idx

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
        outputFilename: bwa_normal_output_filename
        readGroup: bwa_readGroup
        reads: bwa_normal_reads 
        reference: bwa_ref
        threads: bwa_threads

    sam-view-normal:
        tool: samtools-view
        input: bwa-mem-normal/output
        outputName: sam_normal_output_name
        
    picard-sortsam-normal:
        tool: picard-sortsam
        inputFileName_sortSam: sam-view-normal/output
        outputFileName_sortSam: sortsam_normal_outputFileName_sortSam
        tmpdir: sortsam_tmpdir
        validation_stringency: sortsam_validation_stringency
        maxRecordsInRam: sortsam_maxRecordsInRam
        
    picard-markdup-normal:
        tool: picard-markdups
        inputFileName_markDups: picard-sortsam-normal/out
           # source: picard-sortsam-normal/output
           # valueFrom: ${return [ self ];}
        outputFileName_markDups: markdups_normal_outputFileName_markDups
        metricsFile: markdups_normal_metricsFile
        picard_markdup_tmpdir: markdups_picard_markdup_tmpdir
        maxRecordsInRam: markdups_maxRecordsInRam
        validation_stringency: markdups_validation_stringency

    gather-md-normal-index:
        tool: gather
        bamFile: picard-markdup-normal/out
        bamIndex: picard-markdup-normal/out_idx
    
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
        bamFile: gatk-printreads-normal/out
        bamIndex: gatk-printreads-normal/out_idx
        
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
        tumor: gather-pr-tumor-index/output
        normal: gather-pr-normal-index/output
        reference: mutect_reference
        bedFile: mutect_bedFile
        outputFilename: mutect_outputfile_name
        dbsnp: mutect_dbsnp
        cosmic: mutect_cosmic
    
    #==================    
        
        
        
        
outputs:
    o1: bwa-mem-tumor/bwa_output
    o2: bwa-mem-normal/bwa_output
    o3: sam-view-tumor/sam_output
    o4: sam-view-normal/sam_output
    o5: picard-sortsam-tumor/out
    o6: picard-sortsam-tumor/indexes
    o7: picard-sortsam-tumor/out
    o8: picard-sortsam-tumor/indexes
    o9: picard-markdup-tumor/out
    o10: picard-markdup-tumor/metrics
    o11: picard-markdup-tumor/out_idx
    o12: picard-markdup-normal/out
    o13: picard-markdup-normal/metrics
    o14: picard-markdup-normal/out_idx
    o15: gatk-recalibrator-tumor/output_baseRecalibrator
    o16: gatk-recalibrator-normal/output_baseRecalibrator
    o17: gatk-printreads-tumor/out
    o18: gatk-printreads-tumor/out_idx
    o19: gatk-printreads-normal/out
    o20: gatk-printreads-normal/out_idx
    o21: gatk-haplotypecaller-tumor/out
    o22: gatk-haplotypecaller-tumor/bamOut
    o23: gatk-haplotypecaller-normal/out
    o24: gatk-haplotypecaller-normal/bamOut
    o25: gatk-mutect/output_mutect2
    o26: gather-md-normal-index/combined


"""


class ParallelPipeline():

    @staticmethod
    def test_simple():
        from pipeline_definition.pipeline_translator import PipelineTranslator
        Logger.set_write_location("/Users/franklinmichael/source/wehi-pipeline-definition/parallel.log")
        translator = PipelineTranslator("parallel")
        translation = translator.translate_string(_yml)
        Logger.close_file()
        print(translation)


if __name__ == '__main__':
    ParallelPipeline.test_simple()
