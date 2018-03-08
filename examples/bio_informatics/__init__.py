from pipeline_definition.types.type_registry import register_input_factory
from pipeline_definition.types.type_registry import register_step_factory


from examples.bio_informatics.bam import BAMFactory
from examples.bio_informatics.paired_read import PairedReadFactory

from examples.bio_informatics.align import AlignFactory
from examples.bio_informatics.trim import TrimFactory
from examples.bio_informatics.call import CallFactory
from examples.bio_informatics.joint_call import JointCallFactory
from examples.bio_informatics.dedup import DedupFactory

register_input_factory(PairedReadFactory())
register_input_factory(BAMFactory())

register_step_factory(AlignFactory())
register_step_factory(TrimFactory())
register_step_factory(CallFactory())
register_step_factory(JointCallFactory())
register_step_factory(DedupFactory())



