from pipeline_definition.types.type_registry import register_input_factory
from pipeline_definition.types.type_registry import register_step_factory

from examples.bio_informatics.data_types.bam import BamFactory
from examples.bio_informatics.data_types.sorted_bam import SortedBamFactory
from examples.bio_informatics.data_types.reference import ReferenceFactory
from examples.bio_informatics.data_types.paired_read import PairedReadFactory
from examples.bio_informatics.data_types.trimmed_reads import TrimmedReadFactory
from examples.bio_informatics.data_types.text import TextFactory

from examples.bio_informatics.steps.index_bam import IndexBamFactory
from examples.bio_informatics.steps.sort_bam import SortBamFactory
from examples.bio_informatics.steps.align import AlignFactory
from examples.bio_informatics.steps.trim import TrimFactory
from examples.bio_informatics.steps.call import CallFactory
from examples.bio_informatics.steps.joint_call import JointCallFactory
from examples.bio_informatics.steps.dedup import DedupFactory
from examples.bio_informatics.steps.fastqc import FastQCFactory
from examples.bio_informatics.steps.intersect_genic import IntersectFactory

register_input_factory(PairedReadFactory())
register_input_factory(TrimmedReadFactory())
register_input_factory(BamFactory())
register_input_factory(SortedBamFactory())
register_input_factory(ReferenceFactory())
register_input_factory(TextFactory())

register_step_factory(AlignFactory())
register_step_factory(IndexBamFactory())
register_step_factory(SortBamFactory())
register_step_factory(TrimFactory())
register_step_factory(CallFactory())
register_step_factory(JointCallFactory())
register_step_factory(DedupFactory())
register_step_factory(FastQCFactory())
register_step_factory(IntersectFactory())



