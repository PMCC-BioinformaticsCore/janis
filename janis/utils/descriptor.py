from weakref import WeakKeyDictionary
from typing import Any
from janis.utils.logger import Logger


class BaseDescriptor(object):
    def __init__(self, default: Any = None):
        self.default: Any = default
        self.data: WeakKeyDictionary = WeakKeyDictionary()

    def __get__(self, instance, owner):
        # Logger.log(f"Getting value from {instance} and {owner}")
        return self.data.get(instance, self.default)

    def __set__(self, instance, value):
        # Logger.log(f"Setting value {value} to {instance}")
        self.data[instance] = value


class GetOnlyDescriptor(object):
    def __init__(self, default: Any = None):
        self.default: Any = default
        self.data: WeakKeyDictionary = WeakKeyDictionary()

    def __get__(self, instance, owner):
        Logger.log(f"Getting readonly value from {instance} and {owner}")
        return self.data.get(instance, self.default)

    def __set__(self, instance, value):
        if not self.data[instance]:
            return Logger.log("Attempting to write to read-only property")

        Logger.log(f"Setting readonlu value {value} to {instance}")
        self.data[instance] = value
