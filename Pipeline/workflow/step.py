from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput
from Pipeline.translations.cwl.cwl import Cwl
from Pipeline.types.filename import Filename
from Pipeline.utils.logger import Logger

CS = Cwl.WORKFLOW.STEP


class Step:
    def __init__(self, identifier: str, tool: Tool, meta: Dict[str, Any]=None,
                 label: str = None, doc: str = None):
        self._identifier: str = identifier
        self.label = label
        self.doc = doc

        self.__tool: Tool = tool
        self.__meta: Dict[str, Any] = meta

    def id(self):
        return self._identifier

    def __str__(self):
        t = self.__tool.id()
        return f"{self.id()}: {t}"

    def input_value(self, tag: str) -> Optional[Any]:
        return self.__meta[tag] if tag in self.__meta else None

    def set_input_value(self, tag: str, value: str):
        l = self._identifier
        Logger.log(f"Updating '{l}': setting '{tag}' -> '{value}'")
        self.__meta[tag] = value

    def requires(self) -> Dict[str, ToolInput]:
        return self.__tool.inputs_map()

    def provides(self) -> Dict[str, ToolOutput]:
        return self.__tool.outputs_map()

    def tool(self) -> Tool:
        return self.__tool

    def cwl(self):
        d = {
            Cwl.WORKFLOW.STEP.kID: self.id(),
            CS.kRUN: f"tools/{self.tool().id()}.cwl",
            CS.kOUT: [o.tag for o in self.tool().outputs()]
        }

        if self.label:
            d[CS.kLABEL] = self.label

        if self.doc:
            d[CS.kDOC] = self.doc

        return d

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        if item in self.tool().inputs_map():
            return f"{self.id()}/{self.tool().inputs_map()[item].tag}"
        if item in self.tool().outputs_map():
            return f"{self.id()}/{self.tool().outputs_map()[item].tag}"

        tags = ", ".join([f"in.{i.tag}" for i in self.tool().inputs()]
                         + [f"out.{o.tag}" for o in self.tool().outputs()])

        raise AttributeError(f"Tool '{self.tool().id()}' has no identifier '{item}' ({tags})")


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
            CS.kIN: self._prepare_cwl_inputs(),
            ** self.step.cwl()
        }

    def _prepare_cwl_inputs(self):
        ds = []
        ins = self.inputs()
        for k in ins:
            inp = ins[k]
            d = {
                CS.STEP_INPUT.kID: k
            }

            if k in self.connection_map:
                d[CS.STEP_INPUT.kSOURCE] = self.connection_map[k][0]

            elif not inp.input_type.optional:
                    raise Exception(f"Error when building connections for step '{self.id()}', "
                                    f"could not find required connection {k}")

            inp_t = self.inputs()[k].input_type
            if isinstance(inp_t, Filename):
                d[CS.STEP_INPUT.kDEFAULT] = inp_t.generated_filename(self.step.id())
            ds.append(d)
        return ds

    def wdl_map(self) -> List[str]:
        q = []
        for inp in self.step.tool().inputs():
            if inp.tag in self.connection_map:
                q.append(f"{inp.tag} = {self.connection_map[inp.tag][0].replace('/', '.')}")
            elif not inp.optional:
                raise Exception(f"Required option '{inp.tag}' for step '{self.id()}' "
                                f"was not found during conversion to WDL")
        return q