#
# Build the WF schema from registered types
#

from pipeline_definition.types import type_registry


def __input_schema():
    ischema = {}

    for factory in type_registry.get_input_factories():
        ischema[factory.type()] = factory.describe()

    return {
        'inputs': {
            'type': 'dict',
            'keyschema': {'type': 'string'},
            'valueschema': {'schema': ischema}
        }
    }


def __step_scheme():
    ischema = {}

    for factory in type_registry.get_step_factories():
        ischema[factory.type()] = factory.describe()

    return {
        'steps': {
            'type': 'list',
            'keyschema': {'type': 'string'},
            'valueschema': {'schema': ischema}
        }
    }


def schema():
    return dict(__input_schema().items() + __step_scheme().items())


