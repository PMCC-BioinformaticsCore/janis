from pipeline_definition.types.step_type import StepFactory

class CallFactory(StepFactory):
    @classmethod
    def describe(cls):
        return {
            'schema': {
                'caller': {
                    'type': 'string'
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        print(">>>>>>>>>>>>>>>> ", meta )
        return None

    @classmethod
    def type(cls):
        return 'call'

    @classmethod
    def label(cls):
        return 'Caller'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def emit(cls):
        return "Fcatory says: " + cls.__name__

