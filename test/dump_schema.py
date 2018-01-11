#
# Test loading the registry and dumping the schema
#
from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.schema import schema
from pipeline_definition.types.type_registry import register_input_factory


#
class PairedReadFactory(InputFactory):

    @classmethod
    def describe(cls):
        return {
            'forward-pattern':  {'type': 'string'},
            'backward-pattern': {'type': 'string'}
        }

    @classmethod
    def build(cls, yml):
        return None

    @classmethod
    def type(cls):
        return 'SequenceReadArchivePaired'


class BAMFactory(InputFactory):

    @classmethod
    def describe(cls):
        return {
            'path':  {'type': 'string'}
        }

    @classmethod
    def build(cls, yml):
        return None

    @classmethod
    def type(cls):
        return 'BAM'


register_input_factory(PairedReadFactory())
register_input_factory(BAMFactory())

print(schema())
