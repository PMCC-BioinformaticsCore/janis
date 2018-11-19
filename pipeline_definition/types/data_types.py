"""
    Custom Types::

    We want to define some standard interface that allows for the implementation of types.
    Types must be supported as Inputs and Outputs, and must be comparable, ie we can compare
    the output of one step to another step. We'll guarantee that

    We'll also may require separate interfaces for CWL and WDL implementations, we may not be
    able to genericise the interface enough to support an automatic conversion

    We are allowed to require that a type must register itself when loaded.

"""
from abc import ABC


class DataType(ABC):

    def __init__(self, optional=False):
        self.optional = optional

    @staticmethod
    def name():
        raise Exception("Subclass MUST override name field")

    @staticmethod
    def doc():
        return """
        Subclasses should override this class to provide additional information on how to
        correctly provide data to the class, what inputs it may have and what other types
        are compatible
        """

    # The following methods don't need to be overriden, but can be

    def id(self):
        if self.optional:
            return f"Optional<{self.name()}>"
        return self.name()

    def identify(self):
        print(self.id())

    def can_receive_from(self, other) -> bool:
        """
        Can this class receive from $other, likely going to be type(a) == type(b)
        :param other:
        :return:
        """
        if not isinstance(other, type(self)):
            return False
        if self.optional:
            return True
        return self.optional == other.optional

    def input_field_from_input(self, meta):
        return None
