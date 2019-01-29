import janis as p
from bioinformatics.data_types.bampair import BamPair
from bioinformatics.data_types.fasta import FastaWithDict
from bioinformatics.tools.bcftools.view.latest import BcfToolsViewLatest
from bioinformatics.tools.common.splitmultiallele import SplitMultiAllele
from bioinformatics.tools.illumina.strelka.strelka_2_9_9 import Strelka_2_9_9


def strelka_pipeline():
    w = p.Workflow("strelka_pipeline")

    ref = p.Input("reference", FastaWithDict())
    bam = p.Input("bam", BamPair())

    s1_strelka = p.Step("strelka", Strelka_2_9_9())
    s2_bcf = p.Step("bcf", BcfToolsViewLatest())
    s3_split = p.Step("splitmultiallele", SplitMultiAllele())

    w.add_edges([
        (bam, s1_strelka),
        (ref, s1_strelka),
        # (s1_strelka.directory, p.Output("directory")),
        (s1_strelka.script, p.Output("script")),
        (s1_strelka.stats, p.Output("stats")),
        (s1_strelka.variants, p.Output("variants"))
    ])

    # Step2
    w.add_edge(s1_strelka.variants, s2_bcf.file)
    w.add_default_value(s2_bcf.applyFilters, ["PASS"])

    # Step3
    w.add_edges([
        (ref, s3_split.reference),
        (s2_bcf, s3_split.input),
        (s3_split, p.Output("split-out"))
    ])

    return w


w = strelka_pipeline()
w.dump_cwl(to_disk=True, with_docker=True)
