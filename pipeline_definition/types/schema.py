#
# Build the WF document schema from registered types
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
            'valueschema': {
                'schema': ischema
            }
        }
    }


def __step_scheme():
    ischema = {}

    #Inputs can be tagged to indicate processing branch
    ischema['tag'] = {
        'type': 'string',
        'required': False
    }

    ischema['input_scope'] = {
        'type': 'list',
        'required': False
    }

    for factory in type_registry.get_step_factories():
        ischema[factory.type()] = factory.describe()

    return {
        'steps': {
            'type': 'list',
            'schema' : {
                'type' : 'dict',
                'keyschema': {
                    'type': 'string'
                },
                'valueschema': {
                    'schema': ischema
                }
            }
        }
    }


def schema():
    #return dict(__input_schema().items() + __step_scheme().items())

    #https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression
    return {**__input_schema(), **__step_scheme()}

