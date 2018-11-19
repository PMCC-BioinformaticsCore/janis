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
from typing import Dict


class DataType(ABC):

    def __init__(self, optional=False):
        self.optional = optional

    def name(self):
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


###################
# Implementations #
###################

class String(DataType):

    def name(self):
        return "String"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Number(DataType):

    def name(self):
        return "Number"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Boolean(DataType):

    def name(self):
        return "Boolean"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class Array(DataType):
    def __init__(self, t: DataType, optional=False):
        if not isinstance(t, DataType):
            raise Exception(f"Type t ({type(t)}) must be an instance of 'Type'")

        self.__t = t
        super().__init__(optional)

    def name(self):
        return f"Array<{self.__t.id()}>"

    def input_field_from_input(self, meta):
        return next(iter(meta.values()))


class File(DataType):

    def name(self):
        return "File"

    def input_field_from_input(self, meta):
        # WDL: "{workflowName}.label" = meta["path"}
        return {
            "class": "File",
            "path": meta["path"]
        }


class TarFile(File):
    def name(self):
        return "TarFile"


if __name__ == "__main__":
    a = String()
    b = String()
    c = String(optional=True)
    d = File()
    e = Array(String())
    f = Array(String(), optional=True)
    g = Array(String(optional=True), optional=True)
    h = File()
    j = TarFile()

    print("a: " + a.id())
    print("b: " + b.id())
    print("c: " + c.id())
    print("d: " + d.id())
    print("e: " + e.id())
    print("f: " + f.id())
    print("g: " + g.id())
    print("h: " + h.id())
    print("j: " + j.id())
    print("")
    print(f"a => b: {str(b.can_receive_from(a))}\t (same type)")
    print(f"a => c: {str(c.can_receive_from(a))}\t (non-optional => optional)")
    print(f"c => a: {str(a.can_receive_from(c))}\t (optional => non-optional)")
    print(f"d => a: {str(a.can_receive_from(d))}\t (mismatch type)")
    print(f"a => d: {str(d.can_receive_from(a))}\t (mismatch type)")
    print(f"h => j: {str(j.can_receive_from(h))}\t (superclass -> subclass)")
    print(f"j => h: {str(h.can_receive_from(j))}\t (subclass -> superclass)")
