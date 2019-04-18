from typing import Dict, Optional, List

from janis.utils.validators import Validators

from janis.utils.logger import Logger, LogLevel
# import cwlgen.cwlgen as cwl
from janis.types.data_types import DataType, NativeTypes
from janis.workflow.step import ToolOutput
from janis.graph.node import Node, NodeTypes


class Input:

    illegal_keywords = ["input"]

    def __init__(self, identifier: str, data_type: DataType, value=None,
                 label: str=None, doc: str=None, default=None):

        if not Validators.validate_identifier(identifier):
            raise Exception(f"The input identifier '{identifier}' was not validated by '{Validators.identifier_regex}' "
                            f"(must start with letters, and then only contain letters, numbers and an underscore)")

        if identifier in Input.illegal_keywords:
            raise Exception(f"The input identifier '{identifier}' is a reserved keyword "
                            f"({', '.join(Input.illegal_keywords)})")

        Logger.log(f"Creating input '{identifier}' with value: '{value}'")
        self._identifier: str = identifier

        self.label = label
        self.doc = doc

        self.value = value
        self.default = default

        if self.default is not None:
            data_type.optional = True

        self.data_type: DataType = data_type

    def id(self):
        return self._identifier

    def cwl_input(self):
        return self.data_type.cwl_input(self.value)

    def wdl_input(self):
        return self.value


class InputNode(Node):

    def __init__(self, inp: Input):
        Node.__init__(self, NodeTypes.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"": ToolOutput(self.input._identifier, self.input.data_type)}

    def inputs(self):
        return None
