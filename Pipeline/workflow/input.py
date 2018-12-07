from typing import Dict

from Pipeline.types.data_types import DataType
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

    def cwl(self):
        return self.data_type.cwl()

    def input_cwl_yml(self):
        if self.data_type.is_prim:
            return self.meta
        else:
            return {
                "class": self.data_type.primitive(),
                "path": self.data_type.get_value_from_meta(self.meta)
            }


class InputNode(Node):

    def __init__(self, inp: Input):
        super().__init__(NodeTypes.INPUT, inp.id())
        self.input: Input = inp

    def outputs(self) -> Dict[str, ToolOutput]:
        # Program will just grab first value anyway
        return {"output": ToolOutput(self.input.label, self.input.data_type)}

    def inputs(self):
        return None

    @staticmethod
    def input_step_tag_name():
        return 'input'

