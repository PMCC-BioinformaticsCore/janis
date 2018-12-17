from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput
from Pipeline.translations.cwl.cwl import Cwl
from Pipeline.types.filename import Filename
from Pipeline.utils.logger import Logger

CS = Cwl.Workflow.Step


class Step:
    def __init__(self, identifier: str, tool: Tool, meta: Optional[Any]=None,
                 label: str = None, doc: str = None):
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

    def set_input_value(self, tag: str, value: str):
        l = self._identifier
        Logger.log(f"Updating '{l}': setting '{tag}' -> '{value}'")
        self.__meta[tag] = value

    def requires(self) -> Dict[str, ToolInput]:
        return self.tool().inputs_map()

    def provides(self) -> Dict[str, ToolOutput]:
        return self.tool().outputs_map()

    def tool(self) -> Tool:
        return self.__tool

    def cwl(self, is_nested_tool=False):
        run_ref = f"{self.tool().id()}.cwl" if is_nested_tool else f"tools/{self.tool().id()}.cwl"
        d = {
            Cwl.Workflow.Step.kID: self.id(),
            CS.kRUN: run_ref,
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

        raise AttributeError(f"Step '{self.id()}' with tool '{self.tool().id()}' has no identifier '{item}' ({tags})")


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

    def cwl(self, is_nested_tool=False):
        """
        :param is_nested_tool: changes the run reference from tools/toolName.cwl -> toolName.cwl
        :return:
        """
        ins = self.inputs()
        scatterable = []

        dd = {
            **self.step.cwl(is_nested_tool=is_nested_tool),
            CS.kIN: []
        }

        for k in ins:
            inp = ins[k]
            d = {CS.StepInput.kID: k}
            if k not in self.connection_map:
                if inp.input_type.optional:
                    continue
                else:
                    raise Exception(f"Error when building connections for step '{self.id()}', "
                                    f"could not find required connection: '{k}'")

            edge = self.connection_map[k]
            d[CS.StepInput.kSOURCE] = edge.source()
            if edge.scatter:
               scatterable.append(k)
            inp_t = self.inputs()[k].input_type
            if isinstance(inp_t, Filename):
                d[CS.StepInput.kDEFAULT] = inp_t.generated_filename(self.step.id())
            dd[CS.kIN].append(d)

        if len(scatterable) > 0:
            if len(scatterable) > 1:
                Logger.info(f"Discovered more than one scatterable field on step '{self.id()}', "
                            f"deciding scatterMethod to be dot_product")
                dd[CS.kSCATTER_METHOD] = CS.ScatterMethod.kDOT_PRODUCT
            dd[CS.kSCATTER] = scatterable
        return dd


    def wdl_map(self) -> List[str]:
        q = []
        for inp in self.step.tool().inputs():
            if inp.tag in self.connection_map:
                q.append(f"{inp.tag} = {self.connection_map[inp.tag][0].replace('/', '.')}")
            elif not inp.optional:
                raise Exception(f"Required option '{inp.tag}' for step '{self.id()}' "
                                f"was not found during conversion to WDL")
        return q