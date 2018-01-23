from pipeline_definition.types.input_type import InputFactory


class PairedReadFactory(InputFactory):
    @classmethod
    def describe(cls):
        return {
            'schema' : {
                'forward-pattern': {'type': 'string'},
                'backward-pattern': {'type': 'string'}
            }
        }

    @classmethod
    def build(cls, yml):
        return None

    @classmethod
    def type(cls):
        return 'SequenceReadArchivePaired'

    @classmethod
    def label(cls):
        return 'paired read files'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def emit(cls):
        return "Translation by " + cls.__name__
