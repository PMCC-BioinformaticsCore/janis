from abc import ABC, abstractmethod
import re
from typing import List, Dict, Optional, Any

from pipeline_definition.types.data_types import DataType


class ToolInput:
    def __init__(self, tag: str, input_type: DataType, position: int=0, prefix: Optional[str]=None, separate_prefix: bool=True):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        self.tag: str = tag
        self.input_type: DataType = input_type
        self.position: int = position
        self.prefix: Optional[str] = prefix
        self.separate_prefix = separate_prefix

    def cwl(self):
        return {
            "inputBinding": {
                "position": self.position
            },
            **self.input_type.cwl()
        }


class ToolOutput:
    def __init__(self, tag: str, output_type: DataType, glob: Optional[str]=None):
        self.tag = tag
        self.output_type: DataType = output_type
        self.glob = glob

    def cwl(self):
        d = {**self.output_type.cwl()}
        if self.glob is not None:
            d["outputBinding"] = {
                "glob": self.glob
            }
        return d


class ToolArgument:

    expr_pattern = "\$\(.*\)"

    def __init__(self, value: str, prefix: str, position: int=0, separate_value_from_prefix=True):
        self.prefix = prefix
        self.value = value
        self.position = position
        self.is_expression = re.match(self.expr_pattern, self.value) is not None
        self.separate_value_from_prefix = separate_value_from_prefix

    def cwl(self):
        d = {
            "position": self.position,
            "valueFrom": self.value
        }
        if not self.separate_value_from_prefix:
            d["separate"] = self.separate_value_from_prefix

        if self.prefix:
            d["prefix"] = self.prefix

        return d


class Tool(ABC):
    """
    Notes:
        - If you're thinking about secondary files, DON'T!! Consider creating a new DataType.
        - This class is similar to how RABIX COMPOSER creates the tools
        - You can subclass and override whichever fields you'd like, including the INPUTS / OUTPUTS!
        - Take note which options you can provide to the ToolInput and ToolOutput.
    """

    @staticmethod
    @abstractmethod
    def tool():
        raise Exception(f"subclass MUST implement 'tool' method")

    @classmethod
    def full_name(cls):
        if cls.version() is not None:
            return f"{cls.tool()}/{cls.version()}"
        return cls.tool()

    @staticmethod
    def version():
        return None

    @staticmethod
    def docker():
        return "ubuntu:latest"

    @staticmethod
    @abstractmethod
    def base_command():
        raise Exception("Subclass MUST implement the 'base_command' method with str or [str] return types")

    def memory(self) -> Optional[int]:
        return None

    def cpus(self) -> Optional[int]:
        return None

    def arguments(self) -> Optional[List[ToolArgument]]:
        return None

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

    def cwl(self) -> Dict[str, Any]:
        d = {
            "class": "CommandLineTool",
            "cwlVersion": "v1.0",
            "baseCommand": self.base_command(),
            "id": self.tool(),
            "label": self.tool(),
        }

        hints = {}
        if self.docker() is not None:
            hints["DockerRequirement"] = {"dockerPull": self.docker()}

        if hints:
            d["hints"] = hints

        inps = {}
        for tool_input in self.inputs():
            inps[tool_input.tag] = tool_input.cwl()

        if self.inputs():
            d["inputs"] = {t.tag: t.cwl() for t in self.inputs()}

        if self.outputs():
            d["outputs"] = {t.tag: t.cwl() for t in self.outputs()}

        if self.arguments():
            d["arguments"] = [a.cwl() for a in self.arguments()]

        return d