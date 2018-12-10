from typing import Dict

from Pipeline.types.data_types import DataType, NativeTypes
from Pipeline.workflow.step import ToolOutput
from Pipeline.graph.node import Node, NodeTypes


class Input:
    def __init__(self, label: str, data_type: DataType, meta=None):
        self.label: str = label
        self.data_type: DataType = data_type
        # Will have type represented by data_type
        self.meta = meta

    def id(self):
        return self.label


class InputNode(Node):

    def __init__(self, inp: Input):
        super().__init__(NodeTypes.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"output": ToolOutput(self.input.label, self.input.data_type)}

    def inputs(self):
        return None

    def cwl(self):
        return self.input.data_type.cwl()

    def input_cwl_yml(self):
        if self.input.data_type.is_prim:
            return self.input.meta
        else:
            return self.input.data_type.input_field_from_input(self.input.meta)
            # return {
            #     "class": NativeTypes.map_to_cwl(self.input.data_type.primitive()),
            #     "path": self.input.data_type.get_value_from_meta(self.input.meta)
            # }