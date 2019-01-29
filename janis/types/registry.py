from janis.utils.registry import Registry
from janis.types.data_types import DataType

from typing import Type, List, Optional


__types_registry = Registry[Type[DataType]]()


def register_type(data_type: Type[DataType]) -> bool:
    type_id = data_type.name().lower()
    return __types_registry.register(type_id, data_type)


def get_type(type_name: str) -> Optional[Type[DataType]]:
    return __types_registry.get(type_name.lower())


def get_types() -> List[Type[DataType]]:
    return __types_registry.objects()
