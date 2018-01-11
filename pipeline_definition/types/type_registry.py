#
# Registry of all input types
#
from pipeline_definition.utils.registry import Registry

__inputs_registry = Registry()
__steps_registry = Registry()


def register_input_factory(factory):
    __inputs_registry.register_factory(factory)


def get_input_factory(type_name):
    return __inputs_registry.factory*type_name


def get_input_factories():
    return __inputs_registry.factories()


def register_step_factory(factory):
    __steps_registry.register_factory(factory)


def get_step_factory(type_name):
    return __steps_registry.factory*type_name


def get_step_factories():
    return __steps_registry.factories()
