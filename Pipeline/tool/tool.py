from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any
import re

from Pipeline.utils.logger import Logger
from Pipeline.types.data_types import DataType
from Pipeline.translations.cwl.cwl import Cwl

ToolType = str


class ToolTypes:
    Workflow: ToolType = "workflow"
    CommandTool: ToolType = "command-tool"
    ExpressionTool: ToolType = "expression-tool"


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(self, value: str, prefix: Optional[str] = None, position: int = 0, separate_value_from_prefix=True):
        self.prefix = prefix
        self.value = value
        self.position = position
        self.is_expression = re.match(self.expr_pattern, self.value) is not None
        self.separate_value_from_prefix = separate_value_from_prefix

        if self.prefix and not self.separate_value_from_prefix and not self.prefix.endswith("="):
            # I don't really know what this means.
            Logger.warn(f"Argument ({self.prefix} {self.value}) is not separating and did not end with ='")

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

    def wdl(self):
        return (self.prefix if self.prefix is not None else "") \
               + (" " if self.separate_value_from_prefix else "") \
               + (self.value if self.value is not None else "")


class ToolInput(ToolArgument):
    def __init__(self, tag: str, input_type: DataType, position: Optional[int] = 0, prefix: Optional[str] = None,
                 separate_value_from_prefix: bool = True, default: Any = None, doc: Optional[str]=None):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        super().__init__("", prefix, position, separate_value_from_prefix)

        if default is not None:
            input_type.optional = True

        self.tag: str = tag
        self.input_type: DataType = input_type
        self.optional = self.input_type.optional
        self.default = default
        self.doc = doc

    def cwl(self):
        CLTI = Cwl.CommandLineTool.Inputs
        d = {
            CLTI.kID: self.tag,
            **self.input_type.cwl()
        }

        if self.position is not None or self.prefix is not None:
            input_binding = {}
            if self.position:
                input_binding[CLTI.InputBinding.kPOSITION] = self.position
            if self.prefix is not None:
                input_binding[CLTI.InputBinding.kPREFIX] = self.prefix
            d[Cwl.Workflow.Input.kINPUT_BINDING] = input_binding

        if self.default is not None:
            d[CLTI.kDEFAULT] = self.default

        if self.doc is not None:
            d[CLTI.kDOC] = self.doc

        return d


class ToolOutput:
    def __init__(self, tag: str, output_type: DataType, glob: Optional[str] = None, doc: Optional[str]=None):
        self.tag = tag
        self.output_type: DataType = output_type
        self.glob = glob
        self.doc = doc

    def cwl(self):
        d = { **self.output_type.cwl() }
        if self.glob is not None:
            d[Cwl.Workflow.Output.kOUTPUT_BINDING] = {
                "glob": self.glob
            }
        return d


class Tool(ABC):
    """
    One of Workflow, CommandLineTool, ExpressionTool* (* unimplemented)
    """

    @classmethod
    def type(cls) -> ToolType:
        print(cls)
        raise Exception("Must implement type() method")

    @abstractmethod
    def id(self) -> str:
        raise Exception("Must implement id() method")

    @abstractmethod
    def inputs(self) -> List[ToolInput]:
        raise Exception("Must implement inputs() method")

    @abstractmethod
    def outputs(self) -> List[ToolOutput]:
        raise Exception("Must implement outputs() method")

    def inputs_map(self) -> Dict[str, ToolInput]:
        return {inp.tag: inp for inp in self.inputs()}

    def outputs_map(self) -> Dict[str, ToolOutput]:
        return {outp.tag: outp for outp in self.outputs()}

    @abstractmethod
    def cwl(self):
        raise Exception("Must implement cwl() method")

    @staticmethod
    def doc() -> Optional[str]:
        return None

    def help(self):
        import inspect
        path = inspect.getfile(self.__class__)

        ins = sorted(self.inputs(), key=lambda i: i.position)

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (t.prefix + ": ") if t.separate_value_from_prefix else t.prefix
            return f"\t\t{t.tag} ({prefix_with_space}{t.input_type.id()}" \
                f"{('=' + str(t.default)) if t.default is not None else ''}): {'' if t.doc is None else t.doc}"

        output_format = lambda t: f"\t\t{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"

        requiredInputs = "\n".join(input_format(x) for x in ins if not x.optional)
        optionalInputs = "\n".join(input_format(x) for x in ins if x.optional)
        outputs = "\n".join(output_format(o) for o in self.outputs())

        return f"""
Pipeline tool: {path} ({self.id()})
NAME
    {self.id()}

DESCRIPTION
    {self.doc() if self.doc is not None else "No documentation provided"}

INPUTS:
REQUIRED:
{requiredInputs}

OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
    """

