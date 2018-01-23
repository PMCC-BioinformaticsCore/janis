from pipeline_definition.types.type_registry import register_input_factory
from pipeline_definition.types.type_registry import register_step_factory


from pipeline_definition.bio_informatics.bam import BAMFactory
from pipeline_definition.bio_informatics.paired_read import PairedReadFactory

from pipeline_definition.bio_informatics.align import AlignFactory
from pipeline_definition.bio_informatics.trim import TrimFactory
from pipeline_definition.bio_informatics.call import CallFactory

register_input_factory(PairedReadFactory())
register_input_factory(BAMFactory())

register_step_factory(AlignFactory())
register_step_factory(TrimFactory())
register_step_factory(CallFactory())



