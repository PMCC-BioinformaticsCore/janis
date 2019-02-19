from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any, Union
import re

from janis.types import Selector
from janis.types.common_data_types import Array
from janis.utils.logger import Logger
from janis.types.data_types import DataType, NativeTypes
# import cwlgen.cwlgen as cwl
from janis.utils.metadata import Metadata

ToolType = str


class ToolTypes:
    Workflow: ToolType = "workflow"
    CommandTool: ToolType = "command-tool"
    ExpressionTool: ToolType = "expression-tool"


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(self, value: Any, prefix: Optional[str] = None, position: Optional[int] = 0,
                 separate_value_from_prefix=None, doc: Optional[str] = None, shell_quote: bool = None):
        self.prefix: Optional[str] = prefix
        self.value = value
        self.position: Optional[int] = position
        self.is_expression = isinstance(self.value, Selector) \
                             or (re.match(self.expr_pattern, self.value) is not None) if self.value else None
        self.separate_value_from_prefix = separate_value_from_prefix
        self.doc = doc
        self.shell_quote = shell_quote

        if self.prefix and self.separate_value_from_prefix is not None \
                and not self.separate_value_from_prefix and not self.prefix.endswith("="):
            # I don't really know what this means.
            Logger.warn(f"Argument ({self.prefix} {self.value}) is not separating and did not end with ='")


class ToolInput(ToolArgument):
    def __init__(self, tag: str, input_type: DataType, position: Optional[int] = None, prefix: Optional[str] = None,
                 separate_value_from_prefix: bool = None, default: Any = None, doc: Optional[str] = None,
                 nest_input_binding_on_array: bool = True, shell_quote=None):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        super().__init__(value=None, prefix=prefix, position=position,
                         separate_value_from_prefix=separate_value_from_prefix, doc=doc, shell_quote=shell_quote)

        if default is not None:
            input_type.optional = True

        self.tag: str = tag
        self.input_type: DataType = input_type
        self.optional = self.input_type.optional
        self.default = default
        self.nest_input_binding_on_array = nest_input_binding_on_array

    # def cwl(self):
    #     import cwlgen
    #     default = self.default if self.default else self.input_type.default()
    #
    #     data_type = self.input_type.cwl_type()
    #     input_binding = cwlgen.CommandLineBinding(
    #         # load_contents=self.load_contents,
    #         position=self.position,
    #         prefix=self.prefix,
    #         separate=self.separate_value_from_prefix,
    #         # item_separator=self.item_separator,
    #         # value_from=self.value_from,
    #         shell_quote=self.shell_quote,
    #     )
    #
    #     # Binding array inputs onto the console
    #     # https://www.commonwl.org/user_guide/09-array-inputs/
    #     if isinstance(self.input_type, Array) and isinstance(data_type, cwlgen.CommandInputArraySchema):
    #         if self.nest_input_binding_on_array:
    #             input_binding.prefix = None
    #             input_binding.separate = None
    #             nested_binding = cwlgen.CommandLineBinding(
    #                 # load_contents=self.load_contents,
    #                 prefix=self.prefix,
    #                 separate=self.separate_value_from_prefix,
    #                 # item_separator=self.item_separator,
    #                 # value_from=self.value_from,
    #                 shell_quote=self.shell_quote,
    #             )
    #             data_type.inputBinding = nested_binding
    #         else:
    #             input_binding.itemSeparator = ","
    #
    #     return cwlgen.CommandInputParameter(
    #         param_id=self.tag,
    #         label=self.tag,
    #         secondary_files=self.input_type.secondary_files(),
    #         # streamable=None,
    #         doc=self.doc,
    #         input_binding=input_binding,
    #         default=default,
    #         param_type=data_type
    #     )

    def id(self):
        return self.tag


class ToolOutput:
    def __init__(self, tag: str, output_type: DataType, glob: Optional[Union[Selector, str]] = None,
                 doc: Optional[str] = None):
        self.tag = tag
        self.output_type: DataType = output_type
        self.glob = glob
        self.doc = doc

    def id(self):
        return self.tag


class Tool(ABC, object):
    """
    One of Workflow, CommandLineTool, ExpressionTool* (* unimplemented)
    """

    @classmethod
    @abstractmethod
    def type(cls) -> ToolType:
        raise Exception(f"'{cls}' must implement type() method")

    @abstractmethod
    def id(self) -> str:
        raise Exception("Must implement id() method")

    @staticmethod
    def tool_module():
        return None

    @staticmethod
    def tool_provider():
        return None

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
    def friendly_name(self) -> str:
        # maps to CWL label (still exploring for WDL)
        raise Exception("Tools must implement friendly_name() method")

    def metadata(self) -> Optional[Metadata]:
        return None

    @staticmethod
    def version():
        return None

    def doc(self) -> Optional[str]:
        return None

    def help(self):
        import inspect
        path = inspect.getfile(self.__class__)

        ins = sorted(self.inputs(), key=lambda i: (i.position if i.position is not None else 0))

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

        meta = self.metadata() if self.metadata() else Metadata()

        fn = self.friendly_name() if self.friendly_name() else self.id()
        en = f" ({self.id()})" if fn != self.id() else ""

        return f"""
Pipeline tool: {path} ({fn})
NAME
    {fn + en}

DESCRIPTION
    {meta.documentation if meta.documentation else "No documentation provided"}

INPUTS:
REQUIRED:
{requiredInputs}

OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
    """
