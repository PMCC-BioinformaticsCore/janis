from janis import Workflow, Input, String, Step, Array, Logger, Output, File, Filename
from bioinformatics import BamPair
from bioinformatics import Bed
from bioinformatics import Fasta
from bioinformatics import FastaWithDict
from bioinformatics import Fastq
from bioinformatics import VcfIdx
from bioinformatics import BwaMem
from bioinformatics import Gatk3RecalibratorLatest as Gatk3Recalibrator
from bioinformatics import Gatk3HaplotypeCallerLatest as Gatk3HaplotypeCaller
from bioinformatics import Gatk3Mutect2Latest as Gatk3Mutect2
from bioinformatics import Gatk3PrintReadsLatest as Gatk3PrintReads
from bioinformatics import PicardMarkDup
from bioinformatics import PicardSortSam
from bioinformatics import SamToolsViewLatest as SamToolsView


class TestParallel:

    @staticmethod
    def parallel_workflow():
        Logger.mute()
        subworkflow = Workflow("subworkflow")

        sub_inp1_bwa_ref = Input("bwa_ref", Fasta())
        sub_inp1_reads = Input("reads", Fastq())
        sub_inp1_bwa_read_group = Input("bwa_readGroup", String())

        sub_inp3_sortsam_valstrin = Input("sortsam_validation_stringency", String())

        sub_inp5_recal_ref = Input("recalibrator_ref", FastaWithDict())
        sub_inp5_recal_known = Input("recalibrator_known", Array(VcfIdx()))
        sub_inp5_recal_bed = Input("recalibrator_bed", Bed())

        sub_inp6_pr_ref = Input("printread_ref", FastaWithDict())
        sub_inp6_pr_bed = Input("printread_bed", Bed())
        sub_inp6_pr_outname = Input("printread_outputFilename", Filename())

        sub_inp7_haplo_ref = Input("haplo_ref", FastaWithDict())
        sub_inp7_haplo_dbsnp = Input("haplo_dbsnp", VcfIdx())
        sub_inp7_haplo_bedFile = Input("haplo_bed", Bed())
        sub_inp7_haplo_bam_outname = Input("haplo_bam_outname", Filename())
        sub_inp7_haplo_outname = Input("haplo_outname", Filename())

        o1 = Output("paired", BamPair())
        o2 = Output("haplotype_output", File())

        step1_bwa_mem = Step("bwa_mem", BwaMem())
        step2_samtools= Step("samtools", SamToolsView())
        step3_sortsam = Step("sortsam", PicardSortSam())
        step4_markdup = Step("markdup", PicardMarkDup())
        step5_recal = Step("gatk_recal", Gatk3Recalibrator())
        step6_printread = Step("gatk_printread", Gatk3PrintReads())
        step7_haplo = Step("gatk_haplo", Gatk3HaplotypeCaller())

        # Inputs -> Step1
        subworkflow.add_edge(sub_inp1_bwa_ref, step1_bwa_mem)
        subworkflow.add_edge(sub_inp1_reads, step1_bwa_mem)
        subworkflow.add_edge(sub_inp1_bwa_read_group, step1_bwa_mem)

        # Step1 -> Step2
        subworkflow.add_edge(step1_bwa_mem, step2_samtools)

        # Step2 -> Step3
        subworkflow.add_edge(step2_samtools, step3_sortsam)
        subworkflow.add_edge(sub_inp3_sortsam_valstrin, step3_sortsam.validation_stringency)

        # Step3 -> Step4
        subworkflow.add_edge(step3_sortsam.out, step4_markdup)

        # Step4 -> Step5
        subworkflow.add_edge(step4_markdup.outputPair, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_ref, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_known, step5_recal)
        subworkflow.add_edge(sub_inp5_recal_bed, step5_recal)

        # Step5 -> Step6
        subworkflow.add_edge(step5_recal, step6_printread)
        subworkflow.add_edge(step4_markdup.outputPair, step6_printread)
        subworkflow.add_edge(sub_inp6_pr_ref, step6_printread)
        subworkflow.add_edge(sub_inp6_pr_bed, step6_printread)
        subworkflow.add_edge(sub_inp6_pr_outname, step6_printread.outputfile_printReads)

        # Step6 -> Step7
        subworkflow.add_edge(step6_printread.pairedOutput, step7_haplo)
        subworkflow.add_edge(sub_inp7_haplo_ref, step7_haplo)
        subworkflow.add_edge(sub_inp7_haplo_dbsnp, step7_haplo)
        subworkflow.add_edge(sub_inp7_haplo_bedFile, step7_haplo)
        subworkflow.add_edge(sub_inp7_haplo_outname, step7_haplo.outputFilename)
        subworkflow.add_edge(sub_inp7_haplo_bam_outname, step7_haplo.bamOutput)

        # Outputs
        subworkflow.add_edge(step6_printread.pairedOutput, o1)
        subworkflow.add_edge(step7_haplo.output, o2)

        Logger.unmute()
        print(subworkflow.help())

        subworkflow.draw_graph()
        # subworkflow.dump_cwl(to_disk=True)

        return subworkflow

    def test_parallel(self):
        w = Workflow("parallel")

        norm_bwa_reads = Input("normal_bwa_reads", Fastq())
        norm_printread_outputFilename = Input("normal_printread_outputFilename", String())
        norm_haplo_bam_outname = Input("normal_haplo_bam_outname", String())
        norm_haplo_outname = Input("normal_haplo_outname", String())
        norms = [norm_bwa_reads, norm_haplo_bam_outname, norm_haplo_outname, norm_printread_outputFilename]

        tum_bwa_reads = Input("tumor_bwa_reads", Fastq())
        tum_printread_outputFilename = Input("tumor_printread_outputFilename", String())
        tum_haplo_bam_outname = Input("tumor_haplo_bam_outname", String())
        tum_haplo_outname = Input("tumor_haplo_outname", String())
        tums = [tum_bwa_reads, tum_printread_outputFilename, tum_haplo_bam_outname, tum_haplo_outname]

        bwa_readGroup = Input("bwa_readGroup", String())
        bwa_ref = Input("bwa_ref", Fasta())
        sortsam_validation_stringency = Input("sortsam_validation_stringency", String())
        recalibrator_ref = Input("recalibrator_ref", FastaWithDict())
        recalibrator_known = Input("recalibrator_known", Array(VcfIdx()))
        recalibrator_bed = Input("recalibrator_bed", Bed())
        printread_ref = Input("printread_ref", FastaWithDict())
        printread_bed = Input("printread_bed", Bed())
        haplo_ref = Input("haplo_ref", FastaWithDict())
        haplo_dbsnp = Input("haplo_dbsnp", VcfIdx())
        haplo_bed = Input("haplo_bed", Bed())
        mutect_ref = Input("mutect_ref", FastaWithDict())
        mutect_bed = Input("mutect_bed", Bed())
        mutect_dbsnp = Input("mutect_dbsnp", Array(VcfIdx()))
        mutect_cosmic = Input("mutect_cosmic", Array(VcfIdx()))

        sub_normal: Step = Step("normal_subworkflow", self.parallel_workflow())
        sub_tumor: Step = Step("tumor_subworkflow", self.parallel_workflow())
        step_gatk_mutect2 = Step("gatk_mutect", Gatk3Mutect2())

        out_sub_tumor = Output("tumor", BamPair())
        out_sub_normal = Output("normal", BamPair())
        out_sub_mutect = Output("mutect", File())

        # connect dependencies for both subworkflows
        for sub, inps in [[sub_normal, norms], [sub_tumor, tums]]:
            w.add_edges([
                (bwa_readGroup, sub.bwa_readGroup),
                (bwa_ref, sub),     # should work fine I think
                (sortsam_validation_stringency, sub.sortsam_validation_stringency),
                (recalibrator_ref, sub.recalibrator_ref),
                (recalibrator_known, sub.recalibrator_known),
                (recalibrator_bed, sub.recalibrator_bed),
                (printread_ref, sub.printread_ref),
                (printread_bed, sub.printread_bed),
                (haplo_ref, sub.haplo_ref),
                (haplo_dbsnp, sub.haplo_dbsnp),
                (haplo_bed, sub.haplo_bed),

                # specific inputs
                (inps[0], sub.reads),
                (inps[1], sub.printread_outputFilename),
                (inps[2], sub.haplo_bam_outname),
                (inps[3], sub.haplo_outname)
            ])

        w.add_edges([
            (sub_normal, step_gatk_mutect2.inputBam_normal),
            (sub_tumor, step_gatk_mutect2.inputBam_tumor),
            (mutect_ref, step_gatk_mutect2),
            (mutect_bed, step_gatk_mutect2.bedFile),
            (mutect_dbsnp, step_gatk_mutect2.dbsnp),
            (mutect_cosmic, step_gatk_mutect2.cosmic)
        ])

        w.add_edges([
            (sub_normal.paired, out_sub_normal),
            (sub_tumor.paired, out_sub_tumor),
            (step_gatk_mutect2, out_sub_mutect)
        ])

        return w


if __name__ == "__main__":
    w = TestParallel().test_parallel()
    # w.draw_graph()
    w.dump_cwl(to_disk=True, with_docker=True)
