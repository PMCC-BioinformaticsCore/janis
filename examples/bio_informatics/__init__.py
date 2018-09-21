from pipeline_definition.types.type_registry import register_input_factory
from pipeline_definition.types.type_registry import register_step_factory


from examples.bio_informatics.bam import BAMFactory
from examples.bio_informatics.reference import ReferenceFactory
from examples.bio_informatics.paired_read import PairedReadFactory
from examples.bio_informatics.text import TextFactory

from examples.bio_informatics.align import AlignFactory
from examples.bio_informatics.index_bam import IndexBamFactory
from examples.bio_informatics.sort_bam import SortBamFactory
from examples.bio_informatics.align import AlignFactory
from examples.bio_informatics.trim import TrimFactory
from examples.bio_informatics.call import CallFactory
from examples.bio_informatics.joint_call import JointCallFactory
from examples.bio_informatics.dedup import DedupFactory
from examples.bio_informatics.fastqc import FastQCFactory
from examples.bio_informatics.intersect_genic import IntersectFactory

register_input_factory(PairedReadFactory())
register_input_factory(BAMFactory())
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



