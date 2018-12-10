from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput
from Pipeline.types.filename import Filename
from Pipeline.utils.logger import Logger


class Step:
    def __init__(self, label: str, tool: Tool, meta: Dict[str, Any]=None):
        self.__label: str = label
        self.__tool: Tool = tool
        self.__meta: Dict[str, Any] = meta

    def __str__(self):
        l = self.__label
        t = self.__tool
        return f"{l}: {t}"

    def id(self):
        return self.__label

    def input_value(self, tag: str) -> Optional[Any]:
        return self.__meta[tag] if tag in self.__meta else None

    def set_input_value(self, tag: str, value: str):
        l = self.__label
        Logger.log(f"Updating '{l}': setting '{tag}' -> '{value}'")
        self.__meta[tag] = value

    def requires(self) -> Dict[str, ToolInput]:
        return self.__tool.inputs_map()

    def provides(self) -> Dict[str, ToolOutput]:
        return self.__tool.outputs_map()

    def get_tool(self) -> Tool:
        return self.__tool

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        if item in self.get_tool().inputs_map():
            return f"{self.id()}/{self.get_tool().inputs_map()[item].tag}"
        if item in self.get_tool().outputs_map():
            return f"{self.id()}/{self.get_tool().outputs_map()[item].tag}"

        raise AttributeError(f" tool '{self.get_tool().id()}' has no identifier '{item}'")


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
            return f"{self.id()}/{item}"

        raise AttributeError(f"type object '{type(self)}' has no attribute '{item}'")

    def cwl(self):
        return {
            "label": self.step.id(),
            "run": f"tools/{self.step.get_tool().tool().lower()}.cwl",
            "in": {i: self.connection_map[i][0] for i in self.connection_map},
            "out": [o.tag for o in self.step.get_tool().outputs()]
        }

    def wdl_map(self) -> List[str]:
        q = []
        for inp in self.step.get_tool().inputs():
            if inp.tag in self.connection_map:
                q.append(f"{inp.tag} = {self.connection_map[inp.tag][0].replace('/', '.')}")
            elif not inp.optional:
                raise Exception(f"Required option '{inp.tag}' for step '{self.id()}' "
                                f"was not found during conversion to WDL")
        return q


# Should be able to remove the below


class StepFactory(ABC):
    @classmethod
    @abstractmethod
    def type(cls) -> str:
        pass

    @classmethod
    @abstractmethod
    def label(cls) -> str:
        pass

    @classmethod
    def description(cls) -> str:
        return cls.label()

    @classmethod
    @abstractmethod
    def schema(cls) -> dict:
        pass

    @classmethod
    @abstractmethod
    def build(cls, label: str, meta: dict) -> Step:
        pass

    @classmethod
    def support_translations(cls) -> List[str]:
        return ['cwl']

    @classmethod
    def build_from(cls, label: str, step_meta: dict) -> Step:
        step_type = cls.type()
        Logger.log(f"{step_type} factory: Building from {step_meta}")
        obj = cls.build(label, step_meta)
        # obj.identify()
        return obj


class TaggedDatum(ABC):
    @abstractmethod
    def tags(self):
        # A set tags that can select among similar types
        pass

    @abstractmethod
    def datum_type(self):
        # A datum_type
        pass

    def satisfies(self, datum):
        # A concrete implementation here
        pass
