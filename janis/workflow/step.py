from typing import Dict, Any, Optional

from janis.graph.node import Node, NodeTypes
from janis.tool.tool import Tool, ToolInput, ToolOutput


class Step:
    def __init__(
        self,
        identifier: str,
        tool: Tool,
        meta: Optional[Any] = None,
        label: str = None,
        doc: str = None,
    ):
        self._identifier: str = identifier
        self.label = label
        self.doc = doc

        self.__tool: Tool = tool
        self.__meta: Optional[Any] = meta

    def id(self):
        return self._identifier

    def __str__(self):
        t = self.__tool.id()
        return f"{self.id()}: {t}"

    def __repr__(self):
        return str(self)

    def requires(self) -> Dict[str, ToolInput]:
        return self.tool().inputs_map()

    def provides(self) -> Dict[str, ToolOutput]:
        return self.tool().outputs_map()

    def tool(self) -> Tool:
        # haha don't mispell this, otherwise infinite recursion through __getattr__
        return self.__tool

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        if item in self.tool().inputs_map():

            return f"{self.id()}/{self.tool().inputs_map()[item].tag}", self
        if item in self.tool().outputs_map():
            return f"{self.id()}/{self.tool().outputs_map()[item].tag}", self

        tags = ", ".join(
            [f"in.{i.tag}" for i in self.tool().inputs()]
            + [f"out.{o.tag}" for o in self.tool().outputs()]
        )

        raise AttributeError(
            f"Step '{self.id()}' with tool '{self.tool().id()}' has no identifier '{item}' ({tags})"
        )


class StepNode(Node):
    def __init__(self, step: Step):
        super().__init__(NodeTypes.TASK, step.id())
        self.step = step

    def inputs(self) -> Dict[str, ToolInput]:
        return self.step.requires()

    def outputs(self) -> Dict[str, ToolOutput]:
        return self.step.provides()

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        if item in self.inputs() or item in self.outputs():
            return f"{self.id()}/{item}", self

        raise AttributeError(f"type object '{type(self)}' has no attribute '{item}'")

    def __repr__(self):
        return f"{self.node_type}: {self.step.tool().id()}"
