import re
from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any, Union

from janis.types import Selector
from janis.types.data_types import DataType
from janis.utils.logger import Logger
from janis.utils.metadata import Metadata
from janis.utils.validators import Validators

ToolType = str


class ToolTypes:
    Workflow: ToolType = "workflow"
    CommandTool: ToolType = "command-tool"
    ExpressionTool: ToolType = "expression-tool"


class ToolArgument:
    expr_pattern = "\$\(.*\)"

    def __init__(
        self,
        value: Any,
        prefix: Optional[str] = None,
        position: Optional[int] = 0,
        separate_value_from_prefix=None,
        doc: Optional[str] = None,
        shell_quote: bool = None,
    ):
        self.prefix: Optional[str] = prefix
        self.value = value
        self.position: Optional[int] = position
        self.is_expression = (
            isinstance(self.value, Selector)
            or (re.match(self.expr_pattern, self.value) is not None)
            if self.value
            else None
        )
        self.separate_value_from_prefix = separate_value_from_prefix
        self.doc = doc
        self.shell_quote = shell_quote

        if (
            self.prefix
            and self.separate_value_from_prefix is not None
            and not self.separate_value_from_prefix
            and not self.prefix.endswith("=")
        ):
            # I don't really know what this means.
            Logger.warn(
                f"Argument ({self.prefix} {self.value}) is not separating and did not end with ='"
            )


class ToolInput(ToolArgument):
    illegal_keywords = ["input"]

    def __init__(
        self,
        tag: str,
        input_type: DataType,
        position: Optional[int] = None,
        prefix: Optional[str] = None,
        separate_value_from_prefix: bool = None,
        default: Any = None,
        doc: Optional[str] = None,
        prefix_applies_to_all_elements: bool = None,
        shell_quote=None,
        separator=None,
        localise_file=None,
    ):
        """
        :param tag: tag for input, what the yml will reference (eg: input1: path/to/file)
        :param input_type:
        """
        super().__init__(
            value=None,
            prefix=prefix,
            position=position,
            separate_value_from_prefix=separate_value_from_prefix,
            doc=doc,
            shell_quote=shell_quote,
        )

        # if default is not None:
        #     input_type.optional = True

        if not Validators.validate_identifier(tag):
            raise Exception(
                f"The identifier '{tag}' was not validated by '{Validators.identifier_regex}' "
                f"(must startpip instal with letters, and then only contain letters, numbers and an underscore)"
            )

        if tag in ToolInput.illegal_keywords:
            raise Exception(
                f"The input identifier '{tag}' is a reserved keyword "
                f"({', '.join(ToolInput.illegal_keywords)})"
            )

        self.tag: str = tag
        self.input_type: DataType = input_type
        self.default = default
        self.prefix_applies_to_all_elements = prefix_applies_to_all_elements
        self.separator = separator
        self.localise_file = localise_file

    def id(self):
        return self.tag


class ToolOutput:
    illegal_keywords = ["output"]

    def __init__(
        self,
        tag: str,
        output_type: DataType,
        glob: Optional[Union[Selector, str]] = None,
        doc: Optional[str] = None,
    ):
        if not Validators.validate_identifier(tag):
            raise Exception(
                f"The identifier '{tag}' was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )
        if tag in ToolOutput.illegal_keywords:
            raise Exception(
                f"The output identifier '{tag}' is a reserved keyword "
                f"({', '.join(ToolOutput.illegal_keywords)})"
            )
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

    @abstractmethod
    def translate(
        self, translation: str, with_docker=True, with_resource_overrides=False
    ):
        raise Exception("Subclass must provide implementation for 'translate()' method")

    def help(self):
        import inspect

        path = inspect.getfile(self.__class__)

        ins = sorted(
            self.inputs(), key=lambda i: (i.position if i.position is not None else 0)
        )

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (
                    (t.prefix + ": ") if t.separate_value_from_prefix else t.prefix
                )
            return (
                f"\t\t{t.tag} ({prefix_with_space}{t.input_type.id()}"
                f"{('=' + str(t.default)) if t.default is not None else ''}): {'' if t.doc is None else t.doc}"
            )

        output_format = (
            lambda t: f"\t\t{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"
        )

        requiredInputs = "\n".join(
            input_format(x) for x in ins if not x.input_type.optional
        )
        optionalInputs = "\n".join(
            input_format(x) for x in ins if x.input_type.optional
        )
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
