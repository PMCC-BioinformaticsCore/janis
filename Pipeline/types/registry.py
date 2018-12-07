from Pipeline.utils.registry import Registry
from Pipeline.types.data_types import DataType

from typing import Type, List


__types_registry = Registry[DataType]()


def register_type(data_type: Type[DataType]) -> bool:
    type_id = data_type.name().lower()
    return __types_registry.register(type_id, data_type)


def get_type(type_name: str) -> Type[DataType]:
    return __types_registry.get(type_name.lower())


def get_types() -> List[Type[DataType]]:
    return __types_registry.objects()