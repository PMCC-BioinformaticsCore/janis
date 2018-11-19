from abc import ABC, abstractmethod
from typing import List, Dict

from pipeline_definition.types.data_types import DataType


class ToolInput:
    def __init__(self, tag: str, input_type: DataType):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        self.tag = tag
        self.input_type: DataType = input_type


class ToolOutput:
    def __init__(self, tag: str, output_type: DataType):
        self.tag = tag
        self.output_type: DataType = output_type


class Tool(ABC):

    @staticmethod
    @abstractmethod
    def tool():
        raise Exception(f"subclass MUST implement 'tool' method")

    @staticmethod
    def version():
        return None

    @staticmethod
    @abstractmethod
    def supported_translations() -> List[str]:
        raise Exception(f"subclass MUST implement the 'supported_translations' method")

    @classmethod
    def full_name(cls):
        if cls.version() is not None:
            return f"{cls.tool()}/{cls.version()}"
        return cls.tool()

    @abstractmethod
    def inputs(self) -> List[ToolInput]:
        raise Exception(f"{type(self)} MUST implement the 'inputs' method")

    @abstractmethod
    def outputs(self) -> List[ToolOutput]:
        raise Exception(f"{type(self)} MUST implement the 'outputs' method")

    def inputs_map(self) -> Dict[str, ToolInput]:
        return {inp.tag: inp for inp in self.inputs()}

    def outputs_map(self) -> Dict[str, ToolOutput]:
        return {outp.tag: outp for outp in self.outputs()}
