#
# Registry of all types
#
from typing import List

from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.step_type import StepFactory
from pipeline_definition.utils.registry import Registry

__input_factories_registry = Registry[InputFactory]()
__steps_registry = Registry[StepFactory]()


def register_input_factory(factory: InputFactory):
  __input_factories_registry.register(factory.type().type_name(), factory)


def get_input_factory(type_name: str) -> InputFactory:
  return __input_factories_registry.get(type_name)


def get_input_factories() -> List[InputFactory]:
  return __input_factories_registry.objects()


def register_step_factory(factory: StepFactory):
  __steps_registry.register(factory.type(), factory)


def get_step_factory(type_name: str) -> StepFactory:
  return __steps_registry.get(type_name)


def get_step_factories() -> List[StepFactory]:
  return __steps_registry.objects()
