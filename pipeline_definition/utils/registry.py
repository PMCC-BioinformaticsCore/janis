#
# Registry of objects
#
from typing import Dict, List, Generic, TypeVar


class RegistryException(Exception):
  pass


T = TypeVar('T')


class Registry(Generic[T]):
  def __init__(self):
    self.registry: Dict[str, T] = dict()

  def register(self, name: str, obj: T):
    self.registry[name] = obj

  def objects(self) -> List[T]:
    return self.registry.values()

  def get(self, type_name) -> T:
    return self.registry.get(type_name)
