#
# Build the WF document schema from registered types
#

from typing import Dict, Any


class KEYS:
    INPUTS = "inputs"
    STEPS = "steps"
    OUTPUTS = "outputs"



def __input_schema() -> Dict[str, Any]:
    ischema: Dict[str, Any] = {}

    for factory in type_registry.get_input_factories():
        ischema[factory.type().type_name()] = factory.schema()

    return {
        'inputs': {
            'type': 'dict',
            'keyschema': {'type': 'string'},
            'valueschema': {
                'schema': ischema
            }
        }
    }


def __step_scheme() -> Dict[str, Any]:
    ischema: Dict[str, Any] = {
        'tag': {'type': 'string', 'required': False},
        'input_scope': {'type': 'list', 'required': False}
    }

    for factory in type_registry.get_step_factories():
        ischema[factory.type()] = factory.schema()

    return {
        'steps': {
            'type': 'list',
            'schema': {
                'type': 'dict',
                'keyschema': {
                    'type': 'string'
                },
                'valueschema': {
                    'schema': ischema
                }
            }
        }
    }


def workflow_schema() -> Dict[str, Any]:
    # https://stackoverflow.com/questions/38987/how-to-merge-two-dictionaries-in-a-single-expression
    return {**__input_schema(), **__step_scheme()}
