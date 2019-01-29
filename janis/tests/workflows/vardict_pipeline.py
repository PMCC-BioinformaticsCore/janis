from janis import File, Array, String, Float
from bioinformatics.data_types.bampair import BamPair
from bioinformatics.data_types.bed import Bed
from bioinformatics.data_types.fasta import FastaWithDict
from bioinformatics.data_types.vcf import VcfIdx
from bioinformatics.tools.bcftools.annotate.latest import BcfToolsAnnotateLatest
from bioinformatics.tools.common.splitmultiallele import SplitMultiAllele
from bioinformatics.tools.common.vardict import VarDict
from bioinformatics.tools.gatk4.genotypeconcordance.latest import Gatk4GenotypeConcordanceLatest
from bioinformatics.tools.htslib.bgzip.latest import BGZipLatest
from bioinformatics.tools.htslib.tabix.latest import TabixLatest


def create():
    import janis as p
    w = p.Workflow("vardict_pipeline")

    input_bed = p.Input("input_bed", Bed())
    indexed_bam = p.Input("indexed_bam", BamPair())
    reference = p.Input("reference", FastaWithDict())
    header_lines = p.Input("headerLines", File())
    truth_vcf = p.Input("truthVcf", VcfIdx())
    intervals = p.Input("intervals", Array(File()))
    sample_name = p.Input("sampleName", String())
    allele_freq_threshold = p.Input("allelFreqThreshold", Float())

    step1 = p.Step("vardict", VarDict())
    step2 = p.Step("annotate", BcfToolsAnnotateLatest())
    step3 = p.Step("split", SplitMultiAllele())
    step4 = p.Step("zip", BGZipLatest())
    step5 = p.Step("tabix", TabixLatest())
    step6 = p.Step("concord", Gatk4GenotypeConcordanceLatest())

    # Step1
    w.add_edges([
        (input_bed, step1.input),
        (indexed_bam, step1.indexedBam),
        (reference, step1),
        (sample_name, step1.sampleName),
        (sample_name, step1.var2vcfSampleName),
        (allele_freq_threshold, step1.alleleFreqThreshold),
        (allele_freq_threshold, step1.var2vcfAlleleFreqThreshold)

    ])
    w.add_default_value(step1.chromNamesAreNumbers, True)
    w.add_default_value(step1.vcfFormat, True)
    w.add_default_value(step1.chromColumn, 1)
    w.add_default_value(step1.regStartCol, 2)
    w.add_default_value(step1.geneEndCol, 3)


    # Step2
    w.add_edges([
        (header_lines, step2.headerLines),
        (step1, step2)
    ])

    # Step3 - splitmulti
    w.add_edges([
        (reference, step3.reference),
        (step2, step3)
    ])

    # Step4 - BGZip
    w.add_edge(step3, step4)

    # Step5 - Tabix
    w.add_edge(step4, step5)

    # Step6 - genotypeConcordance
    w.add_edges([
        (truth_vcf, step6.truthVCF),
        (step5, step6.callVCF),
        (intervals, step6.intervals)
    ])

    w.add_default_value(step6.treatMissingSitesAsHomeRef, True)


    w.add_edges([
        (step1.output, p.Output("vardicted")),
        (step5.output, p.Output("tabixed")),
        (step6.summaryMetrics, p.Output("summaryMetrics")),
        (step6.detailMetrics, p.Output("detailMetrics")),
        (step6.contingencyMetrics, p.Output("contingencyMetrics"))
    ])


    w.dump_cwl(to_disk=True, with_docker=False)
    w.dump_wdl(to_disk=True, with_docker=False)


VardictPipeline = create

if __name__ == "__main__":
    create()