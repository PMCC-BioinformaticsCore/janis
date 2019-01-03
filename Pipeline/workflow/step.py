from abc import ABC, abstractmethod
from typing import List, Dict, Any, Optional

from Pipeline.graph.node import Node, NodeTypes
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput
import cwlgen.cwlgen as cwl
from Pipeline.types.common_data_types import Filename
from Pipeline.utils.logger import Logger


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

    # def cwl(self, is_nested_tool=False):
    #     run_ref = f"{self.tool().id()}.cwl" if is_nested_tool else f"tools/{self.tool().id()}.cwl"
    #     d = {
    #         Cwl.Workflow.Step.kID: self.id(),
    #         CS.kRUN: run_ref,
    #         CS.kOUT: [o.tag for o in self.tool().outputs()]
    #     }
    #
    #     if self.label:
    #         d[CS.kLABEL] = self.label
    #
    #     if self.doc:
    #         d[CS.kDOC] = self.doc
    #
    #     return d

    def cwl(self, is_nested_tool=False):
        run_ref = ("{tool}.cwl" if is_nested_tool else "tools/{tool}.cwl").format(tool=self.tool().id())
        step = cwl.WorkflowStep(
            step_id=self.id(),
            run=run_ref,
            label=self.label,
            doc=self.doc,
            scatter=None,           # Filled by StepNode
            scatter_method=None     # Filled by StepNode
        )

        step.out = [cwl.WorkflowStepOutput(output_id=o.tag) for o in self.tool().outputs()]

        return step

    def __getattr__(self, item):
        if item in self.__dict__:
            return self.__dict__[item]

        if item in self.tool().inputs_map():

            return f"{self.id()}/{self.tool().inputs_map()[item].tag}", self
        if item in self.tool().outputs_map():
            return f"{self.id()}/{self.tool().outputs_map()[item].tag}", self

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
            return f"{self.id()}/{item}", self

        raise AttributeError(f"type object '{type(self)}' has no attribute '{item}'")

    def cwl(self, is_nested_tool=False):

        step = self.step.cwl(is_nested_tool)

        ins = self.inputs()
        scatterable = []

        for k in ins:
            inp = ins[k]
            if k not in self.connection_map:
                if inp.input_type.optional:
                    continue
                else:
                    raise Exception(f"Error when building connections for step '{self.id()}', "
                                    f"could not find required connection: '{k}'")

            inp_t = self.inputs()[k].input_type
            edge = self.connection_map[k]
            default = edge.default if edge.default else inp_t.default()
            d = cwl.WorkflowStepInput(
                input_id=inp.tag,
                source=edge.source(),
                link_merge=None,        # this will need to change when edges have multiple source_map
                default=default,
                value_from=None
            )
            if edge.scatter:
                scatterable.append(k)

            step.inputs.append(d)

        if len(scatterable) > 0:
            if len(scatterable) > 1:
                Logger.info("Discovered more than one scatterable field on step '{step_id}', "
                            "deciding scatterMethod to be dot_product".format(step_id=self.id()))
                step.scatterMethod = "dot_product"
            step.scatter = scatterable
        return step

    def wdl_map(self) -> List[str]:
        q = []
        for inp in self.step.tool().inputs():
            if inp.tag in self.connection_map:
                q.append(f"{inp.tag} = {self.connection_map[inp.tag][0].replace('/', '.')}")
            elif not inp.optional:
                raise Exception(f"Required option '{inp.tag}' for step '{self.id()}' "
                                f"was not found during conversion to WDL")
        return q
