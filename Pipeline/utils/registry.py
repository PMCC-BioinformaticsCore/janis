#
# Registry of objects
#
from typing import Dict, List, Generic, TypeVar, Optional


class RegistryException(Exception):
    pass


T = TypeVar('T')

class Registry(Generic[T]):
    def __init__(self):
        self.registry: Dict[str, T] = {}

    def register(self, name: str, obj: T):
        if name in self.registry:
            return False
        self.registry[name] = obj
        return True

    def objects(self) -> List[T]:
        return list(self.registry.values())

    def get(self, type_name) -> Optional[T]:
        return self.registry.get(type_name)

    def __contains__(self, item):
        return item in self.registry
