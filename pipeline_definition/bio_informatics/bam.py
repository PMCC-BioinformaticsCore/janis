from pipeline_definition.types.input_type import InputFactory

class BAMFactory(InputFactory):
    @classmethod
    def describe(cls):
        return {
            'schema' : {
                'path': {'type': 'string'}
            }

        }

    @classmethod
    def build(cls, yml):
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
