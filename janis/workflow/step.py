from typing import List, Dict, Any, Optional

# import cwlgen.cwlgen as cwl
from janis.graph.node import Node, NodeTypes
from janis.graph.stepinput import StepInput
from janis.tool.tool import Tool, ToolInput, ToolOutput
from janis.utils.logger import Logger


class Step:
    def __init__(self, identifier: str, tool: Tool, meta: Optional[Any]=None,
                 label: str = None, doc: str = None):
        self._identifier: str = identifier.lower()
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

    def set_input_value(self, tag: str, value: str):
        l = self._identifier
        Logger.log(f"Updating '{l}': setting '{tag}' -> '{value}'")
        self.__meta[tag] = value

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

    def __repr__(self):
        return f"{self.node_type}: {self.step.tool().id()}"

    # def cwl(self, is_nested_tool=False):
    #
    #     step = self.step.cwl(is_nested_tool)
    #
    #     ins = self.inputs()
    #     scatterable = []
    #
    #     for k in ins:
    #         inp = ins[k]
    #         if k not in self.connection_map:
    #             if inp.input_type.optional:
    #                 continue
    #             else:
    #                 raise Exception(f"Error when building connections for step '{self.id()}', "
    #                                 f"could not find required connection: '{k}'")
    #
    #         inp_t = self.inputs()[k].input_type
    #         edge = self.connection_map[k]
    #         default = edge.default if edge.default else inp_t.default()
    #         d = cwl.WorkflowStepInput(
    #             input_id=inp.tag,
    #             source=edge.source(),
    #             link_merge=None,        # this will need to change when edges have multiple source_map
    #             default=default,
    #             value_from=None
    #         )
    #         if edge.has_scatter():
    #             scatterable.append(k)
    #
    #         step.inputs.append(d)
    #
    #     if len(scatterable) > 0:
    #         if len(scatterable) > 1:
    #             Logger.info("Discovered more than one scatterable field on step '{step_id}', "
    #                         "deciding scatterMethod to be dot_product".format(step_id=self.id()))
    #             step.scatterMethod = "dot_product"
    #         step.scatter = scatterable
    #     return step

    def wdl_map(self) -> List[str]:
        q = []
        for inp in self.step.tool().inputs():
            if inp.tag in self.connection_map:
                si: StepInput = self.connection_map[inp.tag]
                q.append("{tag} = {value}".format(tag=si.ftag, value=si.dotted_source()))
            elif not inp.optional:
                raise Exception(f"Required option '{inp.tag}' for step '{self.id()}' "
                                f"was not found during conversion to WDL")
        return q

    def wdl(self, step_identifier: str, step_alias: str):
        import wdlgen.wdlgen as wdl

        ins = self.inputs()

        # One step => One WorkflowCall. We need to traverse the edge list to see if there's a scatter
        # then we can build up the WorkflowCall / ScatterCall
        scatterable = [self.connection_map[k].dotted_source()
                       for k in self.inputs() if k in self.connection_map and self.connection_map[k].has_scatter()]

        # We need to replace the scatterable key(s) with some random variable, eg: for i in iterable:
        ordered_variable_identifiers = ["i", "j", "k", "x", "y", "z", "a", "b", "c", "ii", "jj", "kk", "xx", "yy", "zz"]
        new_to_old_identifier = {k.dotted_source(): k.dotted_source() for k in self.connection_map.values()
                                 if not isinstance(k.dotted_source(), list)}

        # We'll wrap everything in the scatter block later, but let's replace the fields we need to scatter
        # with the new scatter variable (we'll try to guess one based on the fieldname. We might need to eventually
        # pass the workflow inputs to make sure now conflict will arise.
        # Todo: Pass Workflow input tags to wdl scatter generation to ensure scatter var doesn't conflict with inputs
        for s in scatterable:
            new_var = s.split(".")[-1][0].lower()
            while new_var in new_to_old_identifier:
                new_var = ordered_variable_identifiers.pop(0)
            new_to_old_identifier[s] = new_var

        # Let's map the inputs, to the source.
        # We're using a dictionary for the map atm, but WDL requires the format:
        #       fieldName: sourceCall.Output

        inputs_map = {}
        for k in ins:
            inp = ins[k]
            if k not in self.connection_map:
                if inp.input_type.optional:
                    continue
                else:
                    raise Exception(f"Error when building connections for step '{self.id()}', "
                                    f"could not find required connection: '{k}'")

            edge = self.connection_map[k]
            ds = edge.dotted_source()

            if isinstance(ds, list):
                if len(ds) == 1:
                    ds = ds[0]
                elif len(ds) > 1:
                    Logger.critical("Conversion to WDL does not currently support multiple sources")
                    ds = f'[{", ".join(ds)}]'

            if ds in new_to_old_identifier and new_to_old_identifier[ds]:
                inputs_map[k] = new_to_old_identifier[ds]
            elif edge.default is not None:
                if isinstance(edge.default, bool):
                    inputs_map[k] = "true" if edge.default else "false"
                elif isinstance(edge.default, str):
                    inputs_map[k] = f'"{edge.default}"'
                else:
                    inputs_map[k] = edge.default
            else:
                inputs_map[k] = ds

        call = wdl.WorkflowCall(step_identifier, step_alias, inputs_map)

        for s in scatterable:
            call = wdl.WorkflowScatter(new_to_old_identifier[s], s, [call])

        return call

