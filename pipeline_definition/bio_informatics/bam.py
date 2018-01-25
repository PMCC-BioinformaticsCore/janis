from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.input_type import Input

class BAMFactory(InputFactory):
    @classmethod
    def describe(cls):
        return {
            'schema' : {
                'path': {'type': 'string'}
            }

        }

    @classmethod
    def build(cls, meta):
        print(">>>>>>>>>>>>>>>> ", meta )
        return None

    @classmethod
    def type(cls):
        return 'BAM'

    @classmethod
    def label(cls):
        return 'BAM file'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def emit(cls):
        return "Fcatory says: " + cls.__name__


class BAMInput(Input):
    pass
