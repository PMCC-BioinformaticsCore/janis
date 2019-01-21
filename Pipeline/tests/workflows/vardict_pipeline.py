from Pipeline import File, Array
from Pipeline.bioinformatics.data_types.bampair import BamPair
from Pipeline.bioinformatics.data_types.bed import Bed
from Pipeline.bioinformatics.data_types.fasta import FastaWithDict
from Pipeline.bioinformatics.data_types.vcf import VcfIdx
from Pipeline.bioinformatics.tools.bcftools.annotate.latest import BcfToolsAnnotateLatest
from Pipeline.bioinformatics.tools.common.splitmultiallele import SplitMultiAllele
from Pipeline.bioinformatics.tools.common.vardict import VarDict
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.latest import Gatk4GenotypeConcordanceLatest
from Pipeline.bioinformatics.tools.htslib.bgzip.latest import BGZipLatest
from Pipeline.bioinformatics.tools.htslib.tabix.latest import TabixLatest


def create():
    import Pipeline as p
    w = p.Workflow("vardict_pipeline")

    input_bed = p.Input("input_bed", Bed())
    indexed_bam = p.Input("indexed_bam", BamPair())
    reference = p.Input("reference", FastaWithDict())
    header_lines = p.Input("headerLines", File())
    truth_vcf = p.Input("truthVcf", Array(VcfIdx()))
    intervals = p.Input("intervals", Array(File()))

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
        (reference, step1)
    ])
    w.add_default_value(step1.chromNamesAreNumbers, True)
    w.add_default_value(step1.vcfFormat, True)
    w.add_default_value(step1.alleleFreqThreshold, 0.05)
    w.add_default_value(step1.sampleName, "NA12878")
    w.add_default_value(step1.chromColumn, 1)
    w.add_default_value(step1.regStartCol, 2)
    w.add_default_value(step1.geneEndCol, 3)
    w.add_default_value(step1.var2vcfSampleName, "NA12878")
    w.add_default_value(step1.var2vcfAlleleFreqThreshold, 0.05)


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

    w.dump_cwl(to_disk=True, with_docker=False)

VardictPipeline = create

if __name__ == "__main__":
    create()