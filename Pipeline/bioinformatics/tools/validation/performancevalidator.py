from Pipeline import Step, Workflow, Input, Output, Array
from Pipeline.bioinformatics.data_types.vcf import Vcf, VcfIdx
from Pipeline.bioinformatics.tools.gatk4.genotypeconcordance.latest import Gatk4GenotypeConcordanceLatest
from Pipeline.bioinformatics.tools.htslib.bgzip.bgzip_1_2_1 import BGZip_1_2_1
from Pipeline.bioinformatics.tools.htslib.tabix.tabix_1_2_1 import Tabix_1_2_1


def create_performance_validator_1_2_1():

    w: Workflow = Workflow("performanceValidator")

    inp = Input("vcf", Vcf())
    inp_truth = Input("truth", VcfIdx())
    inp_intervals = Input("intervals", Array(VcfIdx()))

    s1_bgzip = Step("s1_bgzip", BGZip_1_2_1())
    s2_tabix = Step("s2_tabix", Tabix_1_2_1())
    s3_genotypeconcord = Step("s3_genotypeconcord", Gatk4GenotypeConcordanceLatest())

    w.add_pipe(inp, s1_bgzip, s2_tabix, s3_genotypeconcord)
    w.add_edges([(inp_truth, s3_genotypeconcord.truthVCF), (inp_intervals, s3_genotypeconcord.intervals)])
    w.add_default_value(s3_genotypeconcord.treatMissingSitesAsHomeRef, True)

    w.add_edges([
        (s3_genotypeconcord.summaryMetrics, Output("summaryMetrics")),
        (s3_genotypeconcord.detailMetrics, Output("detailMetrics")),
        (s3_genotypeconcord.contingencyMetrics, Output("contingencyMetrics"))
    ])

    return w


PerformanceValidator_1_2_1 = create_performance_validator_1_2_1



if __name__ == "__main__":
    print(PerformanceValidator_1_2_1().help())