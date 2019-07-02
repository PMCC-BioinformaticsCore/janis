from typing import Dict

from janis.graph.node import Node, NodeTypes
from janis.types.data_types import DataType
from janis.utils.logger import Logger
from janis.utils.validators import Validators
from janis.workflow.step import ToolOutput


class Input:

    illegal_keywords = ["input"]

    def __init__(
        self,
        identifier: str,
        data_type: DataType,
        value=None,
        label: str = None,
        doc: str = None,
        default=None,
        include_in_inputs_file_if_none=True,
    ):

        if not Validators.validate_identifier(identifier):
            raise Exception(
                f"The input identifier '{identifier}' was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )

        if identifier in Input.illegal_keywords:
            raise Exception(
                f"The input identifier '{identifier}' is a reserved keyword "
                f"({', '.join(Input.illegal_keywords)})"
            )

        Logger.log(f"Creating input '{identifier}' with value: '{value}'")
        self._identifier: str = identifier

        self.label = label
        self.doc = doc

        self.default = default
        self.include_in_inputs_file_if_none = include_in_inputs_file_if_none

        if self.default is not None:
            data_type.optional = True

        self.data_type: DataType = data_type

        self.value = value
        if not self.validate_value(True):
            raise TypeError(
                f"Value '{str(value)}' for Input '{identifier}' was not valid for '{data_type.id()}' type"
            )

    def id(self):
        return self._identifier

    def cwl_input(self):
        return self.data_type.cwl_input(self.value)

    def wdl_input(self):
        return self.value

    def validate_value(self, allow_null_if_not_optional: bool) -> bool:
        return self.default or self.data_type.validate_value(
            self.value, allow_null_if_not_optional
        )


class InputNode(Node):
    def __init__(self, inp: Input):
        Node.__init__(self, NodeTypes.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"": ToolOutput(self.input._identifier, self.input.data_type)}

    def inputs(self):
        return None
