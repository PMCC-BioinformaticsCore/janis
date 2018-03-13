from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.types.step_type import Step

class IntersectFactory(StepFactory):
    @classmethod
    def type(cls):
        return 'intersect'

    @classmethod
    def label(cls):
        return 'intersect'

    @classmethod
    def description(cls):
        return cls.label()

    @classmethod
    def describe(cls):
        return {
            'schema': {
                "split": {
                    "type": "boolean",
                    "default": false
                },
                "reportNoOverlaps": {
                    "type": "boolean",
                    "default": false
                }
            },
            'nullable': True
        }

    @classmethod
    def build(cls, meta):
        step = IntersectStep( meta )
        return step

class IntersectStep(Step):

    def provides(self):
        return [
            {
                Step.STR_ID : "reports",
                Step.STR_TYPE: "Text"
            }
        ]



    def requires(self):
        return [
            {
                Step.STR_ID: "read",
                Step.STR_TYPE: "SequenceReadArchivePaired"
            }
        ]

