#
# Registry of all types
#
from typing import List, Type

# from pipeline_definition.types.input_type import InputFactory
from pipeline_definition.types.step import Tool
from pipeline_definition.types.data_types import DataType
from pipeline_definition.utils.registry import Registry

# __input_factories_registry = Registry[InputFactory]()
# __steps_registry = Registry[StepFactory]()
__tools_registry = Registry[Tool]()
__types_registry = Registry[DataType]()


# INPUT FACTORIES

# def register_input_factory(factory: InputFactory):
#     __input_factories_registry.register(factory.type().type_name(), factory)
#
#
# def get_input_factory(type_name: str) -> InputFactory:
#     return __input_factories_registry.get(type_name)
#
#
# def get_input_factories() -> List[InputFactory]:
#     return __input_factories_registry.objects()


# STEP FACTORIES

# def register_step_factory(factory: StepFactory):
#     __steps_registry.register(factory.type(), factory)
#
#
# def get_step_factory(type_name: str) -> StepFactory:
#     return __steps_registry.get(type_name)
#
#
# def get_step_factories() -> List[StepFactory]:
#     return __steps_registry.objects()


# TOOLS

def register_tool(tool: Type[Tool]):
    __tools_registry.register(tool.tool().lower(), tool)


def get_tool(tool_name: str) -> Type[Tool]:
    return __tools_registry.get(tool_name.lower())


def get_tools() -> List[Type[Tool]]:
    return __tools_registry.objects()

# Types


def register_type(data_type: Type[DataType]):
    __types_registry.register(data_type.name().lower(), data_type)


def get_type(type_name: str) -> Type[DataType]:
    return __types_registry.get(type_name.lower())


def get_types() -> List[Type[DataType]]:
    return __types_registry.objects()
