import io
import yaml
from typing import Dict, List, Optional, Any

import Pipeline as pp
from wehi.schema import KEYS

Logger = pp.Logger


class Wehi:

    def __init__(self, name: str):
        self.workflow = pp.Workflow(name)
        self.doc: Dict[str, Any] = {}
        self.inputs: List[pp.Input] = []
        self.steps: List[pp.Step] = []
        self.outputs: List[pp.Output] = []

    def parse_string(self, yml: str):
        f = io.StringIO(yml)
        self.doc = yaml.load(f)
        inputs: Dict[str, Any] = self.doc[KEYS.INPUTS]
        steps: Dict[str, Any] = self.doc[KEYS.STEPS]
        outputs: Dict[str, Any] = self.doc[KEYS.OUTPUTS]

        self.inputs = self.parse_inputs(inputs)
        self.steps = self.parse_steps(steps)
        self.outputs = self.parse_outputs(outputs)
        self.build_graph()
        Logger.log("Built graph")
        # self.workflow.draw_graph()

    def build_graph(self):
        for inp in self.inputs:
            self.workflow.add_input(inp)
        for step in self.steps:
            self.workflow.add_step(step)
        for out in self.outputs:
            self.workflow.add_output(out)

        # Now we'll connect edges
        for step in self.steps:
            for tool_tag in step.tool().inputs():
                inp_tag = step.input_value(tool_tag.tag)
                if not inp_tag:
                    if tool_tag.optional: continue
                    #   2. (b)
                    raise Exception(f"Step '{step.id()}' (tool: '{step.tool().id()}') did not contain"
                                    f" the required input '{tool_tag.tag}' with type: '{tool_tag.input_type.id()}'")

                self.workflow.add_edge(inp_tag, f"{step.id()}/{tool_tag.tag}")

        for out in self.outputs:
            self.workflow.add_edge(out.meta, out)

        print(self.workflow)

    @staticmethod
    def parse_inputs(inputs: Dict[str, Any]) -> List[pp.Input]:
        return [Wehi.parse_input(inp_id, meta) for inp_id, meta in inputs.items()]

    @staticmethod
    def parse_steps(steps: Dict[str, Any]) -> List[pp.Step]:
        return [Wehi.parse_step(step_id, meta) for step_id, meta in steps.items()]

    @staticmethod
    def parse_outputs(outputs: Dict[str, Any]) -> List[pp.Output]:
        return [Wehi.parse_output(out_id, meta) for out_id, meta in outputs.items()]

    @staticmethod
    def parse_input(inp_id: str, meta: Dict[str, Any]) -> pp.Input:
        Logger.log(f"Parsing input: '{inp_id}'")
        input_type: pp.DataType = Wehi._parse_known_type(inp_id, meta)
        Logger.log(f"Detected '{inp_id}' as type: '{input_type.id()}'")

        inp = input_type.get_value_from_meta(meta)
        if inp is None:
            raise Exception(f"Could not find input value for '{inp_id}'")

        return pp.Input(inp_id, input_type, inp)

    @staticmethod
    def parse_step(step_id: str, meta: Dict[str, Any]) -> pp.Step:
        Logger.log(f"Parsing step with id: '{step_id}'")
        step_type: Optional[str]

        if isinstance(meta, str):  # I don't think we'll allow this to happen, because it still
            step_type = str(meta)  # implies type inferencing to some point, or maybe all outputs go to this
            meta = {step_type: None}
        elif isinstance(meta, dict):
            step_type = meta[KEYS.Outputs.TOOL]
        else:
            step_type = None

        if step_type is None:
            raise ValueError(f"There was no tool specified for step: {step_id}")

        tool_type = pp.get_tool(step_type)
        if tool_type is None:
            raise ValueError(f"Could not find the tool '{step_type}' for step: {step_id}")

        # step_obj = step_factory.build_from(step_id, meta)
        step = pp.Step(step_id, tool_type(), meta)

        Logger.log(f"Detected '{step.id()}' with tool '{step.tool().id()}'")
        return step

    @staticmethod
    def parse_output(output_id: str, meta: str) -> pp.Output:
        return pp.Output(output_id, None, meta)

    @staticmethod
    def _parse_known_type(input_id, meta) -> pp.DataType:
        if isinstance(meta, str):
            # tool = pp.get_tool(meta)
            # if tool:
            #     return tool
            return pp.String()
        elif isinstance(meta, int):
            return pp.Int()
        elif isinstance(meta, float):
            return pp.Float()
        elif isinstance(meta, bool):
            return pp.Boolean()
        elif isinstance(meta, list):
            t = Wehi._parse_array_type(meta, input_id)
            return pp.Array(t)
        elif isinstance(meta, dict):
            if "type" not in meta:
                raise Exception(f"Must include 'type' field for input: '{input_id}'")
            detected_type = meta["type"]
            matched_type = pp.get_type(detected_type)
            if matched_type is None:
                raise Exception(f"No input type '{detected_type}' registered for input '{input_id}'")
            return matched_type()
        else:
            raise ValueError(f"Unrecognised type '{type(meta)}' for input {input_id}")

    @staticmethod
    def _parse_array_type(meta: List, input_id) -> Optional[pp.DataType]:
        types = {}
        for x in meta:
            if "type" in x:
                t = pp.get_type(x["type"])()
            else:
                t = Wehi._parse_known_type(input_id, x)
            types[t.name()] = t

        if len(types) == 0:
            Logger.log(f"Could not determine type of array for input: '{input_id}'", pp.LogLevel.WARNING)
            raise Exception(f"Could not determine type for input: '{input_id}'")
        elif len(types) == 1:
            return next(iter(types.values()))

        raise Exception(f"This version does not support multi-typed arrays for input '{input_id}', detected types: "
                        f"{','.join(iter(types.keys()))}. It also does not support determining the "
                        f"lowest common denominator type. Please submit a bug if you think this is in error.")