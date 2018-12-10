"""
    Pipeline Translator

        This class manages the DAG (directed-acyclic graph) for the workflow,
        including creation by the WEHI Pipeline Definition, and some work for
        beginning the conversion to CWL.

    History:
        - 2018-03-07 - Created by Mohammad
        - 2018-09-06 - Initial refactor by Evan
        - 2018-11-08 - Begin typing + refactor by Michael
"""

import os
import io
import json
import yaml
from typing import List, Dict, Any, Optional
from timeit import default_timer as timer

import networkx as nx
# from networkx.drawing.nx_agraph import  import graphviz
from networkx.readwrite import json_graph
from cerberus import Validator

import pipeline_definition.types.schema as wehi_schema
import utils.errors as errors
from Pipeline import String, Number, Boolean, Array
from types.data_types import DataType
from workflow.output import Output, OutputNode
from Tool.tool import ToolInput, Tool

from utils.logger import Logger, LogLevel

from graph.node import Node, NodeType, layout_nodes2
from workflow.input import InputNode
from utils import str_presenter
from pipeline_definition.types.type_registry import get_tool, get_type
# from pipeline_definition.types.type_registry import get_input_factory
from Workflow.step import Step, ToolOutput, StepNode
from workflow.input import Input  # , InputFactory, InputType

yaml.add_representer(str, str_presenter)


# class MappedInput:
#     def __init__(self, inputs: List[InputType], candidates, step_input: InputType):
#         self.inputs = inputs
#         self.candidates = candidates
#         self.input_type = step_input.type_name()
#
#     def __repr__(self):
#         return f'inputs={self.inputs}: candidates={self.candidates}'


class PipelineTranslator:
    def __init__(self, label: str):
        Logger.log("Creating PipelineTranslator", LogLevel.INFO)

        self.label: str = label
        self.__input_step: InputNode = None
        self.__work_graph: nx.Graph = None
        self.__labels_map: Dict[str, Node] = {}

        self.inputs: List[Input] = []
        self.steps: List[Step] = []
        self.outputs: List[Output] = []

    def _dump_as_yaml(self, doc: Dict[str, Any], indent: int = 2, prefix="YAML DOC") -> None:
        """
        Dump the provided dictionary as a YAML document
        :param doc: Dict[str, Any]
        """
        # Diagnostic - what have we got?
        Logger.log(f"{prefix} [")
        Logger.log(yaml.dump(doc, indent=indent))
        Logger.log(f"] END {prefix}")

    def __do_translate(self, doc: Dict[str, Any]) -> None:
        """
        Convert dictionary (of WEHI pipeline definition) into DAG
        :param doc: Dict[str, str]
        :return: N/A
        """
        # Do schema validation
        # self.validate_schema(doc)

        # Now convert the YAML doc?
        self._translate_yaml_doc(doc)
        Logger.log("Translation Output: [")
        Logger.log(self.__work_graph)
        Logger.log("]")

    def validate_schema(self, yaml_doc: Dict[str, Any]) -> bool:
        """
        Validate the yaml_doc against the WEHI spec,
        throws exception if cerberus gives errors.
        :param yaml_doc:
        :return:
        """
        self._dump_as_yaml(yaml_doc)

        sch = wehi_schema.workflow_schema()
        self._dump_as_yaml(sch, indent=2, prefix="PDX SCHEMA")

        v: Validator = Validator(sch)

        validation_success: bool = v.validate(yaml_doc)

        if validation_success:
            Logger.log("Schema is Valid")
        else:
            msg = "ERROR! Pipeline definition document validation failed."
            Logger.log(msg)
            if v.errors is not None:
                raise ValueError(msg, v.errors)
        return validation_success

    def _translate_yaml_doc(self, doc: Dict[str, Any]) -> None:
        """
        Pulls out the inputs, outputs, steps and calls to create workflow_graph
        :param doc: Dict[str, Any]
        """
        # Create a in memory instances for all the inputs, steps and outputs

        Logger.log("Starting translation from converted dictionary", LogLevel.INFO)

        inputs: Dict[str, Any] = doc[wehi_schema.KEYS.INPUTS]
        steps: Dict[str, Any] = doc[wehi_schema.KEYS.STEPS]
        outputs: Dict[str, Any] = doc[wehi_schema.KEYS.OUTPUTS]

        if inputs is None:
            # This is probably okay, depends on whether the the steps reference anything
            # Probably better to be a warning which we can log at the end
            raise errors.InvalidInputsException("There were no inputs provided to the converter")

        if steps is None or not steps:
            raise errors.InvalidStepsException("There are no steps in the diagram")

        # if outputs is None:
        #    raise ValueError("No output?")

        self.inputs: List[Input] = self._build_inputs(inputs)
        self.steps: List[Step] = self._build_steps(steps)
        self.outputs: List[Output] = self._build_outputs(outputs, self.inputs, self.steps)

        self.__work_graph = self._create_workflow_graph(self.steps, self.inputs,
                                                        self.outputs)
        # self.draw_graph()
        self._dump_graph()

        # # Now translate the workflow steps

    def draw_graph(self):
        import matplotlib.pyplot as plt

        default_color = 'black'

        G = self.__work_graph
        edges_attributes = [G.edges[e] for e in G.edges]
        edge_colors = [x["color"] if 'color' in x else default_color for x in edges_attributes]
        node_colors = [NodeType.to_col(x.node_type) for x in G.nodes]

        pos = layout_nodes2(list(G.nodes))

        for n in pos:
            G.node[n]['pos'] = pos[n]

        nx.draw(G, pos=pos, edge_color=edge_colors, node_color=node_colors, with_labels=True)
        plt.show()

    def _build_inputs(self, inputs: Dict[str, Any]) -> List[Input]:
        """
        Build a list of inputs from the inputs dictionary
        :param inputs: Dict[str, Any], where Any = str | { $Type: InputData }
        :return: List[Input] >> All the inputs in the workflow
        """
        input_set: List[Input] = []
        Logger.log("Building inputs", LogLevel.DEBUG)

        for input_id, meta in inputs.items():  # StepLabel, StepProperties/Data

            Logger.log(f"Processing input: '{input_id}'")
            input_type: DataType = self._parse_known_type(meta, input_id)
            Logger.log(f"Detected '{input_id}' as type: '{input_type.id()}'")

            # Build an Input from this data
            # input_obj = inp_factory.build_from(input_id, meta)
            input_obj = Input(input_id, input_type, meta)
            input_set.append(input_obj)

        Logger.log(f"Successfully constructed {len(input_set)} inputs", LogLevel.INFO)

        return input_set

    def _parse_known_type(self, meta, input_id) -> DataType:
        if isinstance(meta, str):
            return String()
        elif isinstance(meta, int) or isinstance(meta, float):
            return Number()
        elif isinstance(meta, bool):
            return Boolean()
        elif isinstance(meta, list):
            # Todo: Be recursive here
            t = self._parse_array_type(meta, input_id)
            return Array(t)
        elif isinstance(meta, dict):
            if "type" not in meta:
                raise Exception(f"Must include 'type' field for input: '{input_id}'")
            detected_type = meta["type"]
            matched_type = get_type(detected_type)
            if matched_type is None:
                raise Exception(f"No input type '{detected_type}' registered for input '{input_id}'")
            return matched_type()
        else:
            raise ValueError(f"Unrecognised type '{type(meta)}' for input {input_id}")

    def _parse_array_type(self, meta: List, input_id) -> Optional[DataType]:
        types = {}
        for x in meta:
            if "type" in x:
                t = get_type(x["type"])()
            else:
                t = self._parse_known_type(x)
            types[t.name()] = t

        if len(types) == 0:
            Logger.log(f"Could not determine type of array for input: '{input_id}'", LogLevel.WARNING)
            return None
        elif len(types) == 1:
            return next(iter(types.values()))

        raise Exception(f"This version does not support multi-typed arrays for input '{input_id}', detected types: "
                        f"{','.join(iter(types.keys()))}. It also does not support determining the "
                        f"lowest common denominator type. Please submit a bug if you think this is in error.")

    @staticmethod
    def _build_outputs(outputs: Dict[str, Any], inputs: List[Input], steps: List[Step]) -> List[Output]:
        """
        We require the inputs / outputs because we need the types
        :param outputs:
        :param inputs:
        :param steps:
        :return:
        """
        outputs_ar: List[Output] = []
        input_map: Dict[str, Input] = {s.id(): s for s in inputs}
        step_map: Dict[str, Step] = {s.id(): s for s in steps}

        for output_id in outputs:

            output_source = outputs[output_id]
            inp_tag_parts = output_source.split("/")
            req_node = inp_tag_parts[0]

            if req_node in input_map:
                raise Exception(f"Not allowed to attach the input '{req_node}' to the output '{output_id}'")
            if req_node not in step_map:
                raise Exception(f"Could not find step '{req_node}' to attach to output '{output_id}")

            # Beauty
            step = step_map[req_node]
            provides: Dict[str, ToolOutput] = step.provides()

            if len(provides) == 0:
                raise Exception(
                    f"The step '{req_node}' does not contain any outputs to connect to the output: '{output_id}'")

            if len(provides) == 1:
                so = provides[next(iter(provides))]
                if (len(inp_tag_parts)) == 1:
                    new_output_source = f"{output_source}/{so.tag}"
                    Logger.log(f"Output '{output_id}' should fully specify step/outputname: "
                               f"{output_source} → {new_output_source}", LogLevel.WARNING)
                    output_source = new_output_source
                    Logger.log(f"Updated output '{output_id}' -> {new_output_source}")
                elif inp_tag_parts[1] not in provides:
                    new_output_source = f"{inp_tag_parts[0]}/{so.tag}"
                    Logger.log(f"Output '{output_id}' did not correctly specify an output of '{req_node}', "
                               f"this has been corrected as the step '{req_node}' only contained one ouput: "
                               f"{output_source} → {new_output_source}", LogLevel.WARNING)
                    Logger.log(f"Updated output '{output_id}' -> {new_output_source}")
                    output_source = new_output_source

                output = Output(output_id, output_source, so.output_type)
                outputs_ar.append(output)

            else:
                if len(inp_tag_parts) == 1:
                    raise Exception(f"Output '{output_id}' ({output_source}) did not fully specify step/outputname")
                if inp_tag_parts[1] not in provides:
                    raise Exception(
                        f"The output '{output_id}' could not find an output called '{inp_tag_parts[1]}' in '{req_node}'")
                so = provides[inp_tag_parts[1]]
                output = Output(output_id, output_source, so.output_type)
                outputs_ar.append(output)
        return outputs_ar

    def _build_steps(self, steps: Dict[str, Any]) -> List[Step]:
        """
        Build list of steps (unconnected) purely from the workflow
        :param steps:
        :return:
        """
        #       For this, we'll need to define a specific format (probably close to CWL)
        #       like: - { tool: type/version, ...inputs }

        pipeline_steps: List[Step] = []

        for step_id, meta in steps.items():

            step_type: Optional[str]

            if isinstance(meta, str):  # I don't think we'll allow this to happen, because it still
                step_type = str(meta)  # implies type inferencing to some point, or maybe all outputs go to this
                meta = {step_type: None}
            elif isinstance(meta, dict):
                step_type = Step.select_type_name_from(meta)
            else:
                step_type = None

            if step_type is None:
                raise ValueError(f"There was no tool specified for step: {step_id}")

            Logger.log(f"Processing STEP: {step_id} with tool {step_type}")

            tool_type = get_tool(step_type)
            if tool_type is None:
                raise ValueError(f"Could not find the tool '{step_type}' for step: {step_id}")

            # step_obj = step_factory.build_from(step_id, meta)
            step_obj = Step(step_id, tool_type(), meta)
            pipeline_steps.append(step_obj)

        return pipeline_steps

    def _create_workflow_graph(self, pipeline_steps: List[Step], workflow_inputs: List[Input],
                               workflow_outputs: List[Output]) -> nx.Graph:
        """
        Take the inputs, {outputs} and steps and stitch it all together
        :param pipeline_steps: List[Step]
        :param workflow_inputs: List[Input]
        :return: a networx graph
        """
        start = timer()
        Logger.log("Building DAG...", LogLevel.INFO)
        work_graph = nx.MultiDiGraph()  # Use the DAG from networx

        # Lets create the input cur_step:
        #   the start cur_step of every workflow that produces the workflow inputs as cur_step output
        # self.__input_step = InputNode(workflow_inputs)
        # pipeline_steps.insert(0, self.__input_step)

        input_nodes = [InputNode(inp) for inp in workflow_inputs]
        work_graph.add_nodes_from(input_nodes)

        step_nodes = [StepNode(step) for step in pipeline_steps]
        work_graph.add_nodes_from(step_nodes)

        output_nodes = [OutputNode(output) for output in workflow_outputs]
        work_graph.add_nodes_from(output_nodes)

        labels = {n.label: n for n in work_graph.nodes}

        # Todo: Move the matching and type checking logic to the Input creation

        for step_node in step_nodes:
            cur_step = step_node.step
            step_tool = cur_step.get_tool()  # Tool that the cur_step references
            required_inputs = cur_step.requires()  # All inputs required for the tool

            for inp_tag in required_inputs:

                # ORDER:
                #   1. Get the input type to link
                #   2. Try and find this type in the cur_step body
                #       (a) - If it's optional, skip
                #       (b) - else: throw error
                #   3. Determine where the input is coming from (input / other cur_step)
                #   4. Try to exactly connect the input to the output/tag
                #       (a) If there's only one output, log warning if it doesn't exactly match
                #       (b) If there are multiple outputs, raise Exception if it doesn't exactly match
                #   5. Check that the types match
                #   6. Add edge
                #

                #   1.
                input_tool: ToolInput = required_inputs[inp_tag]  # The current input we're testing
                input_type: DataType = input_tool.input_type  # The DataType we need

                #   2.
                input_tag: str = cur_step.input_value(inp_tag)  # The tag

                if input_tag is None:
                    #   2. (a)
                    if input_type.optional: continue
                    #   2. (b)
                    raise Exception(f"Step '{step_node.label}' (tool: '{step_tool.tool()}') did not contain"
                                    f" the required input '{inp_tag}' with type: '{input_type.id()}'")

                #   3.
                # Get the correct label, to build the acylic graph, we really only need to first section when split by /
                if type(input_tag) != str:
                    raise Exception(f"Unexpected type {type(input_tag)} for cur_step: '{cur_step.id()}/{inp_tag}'")
                inp_tag_parts = input_tag.split("/")
                required_step_label = inp_tag_parts[0]

                if required_step_label not in labels:
                    raise Exception(f"Could not find node '{required_step_label} when building cur_step '{step_node}'")

                input_node: Node = labels[required_step_label]
                input_node.set_depth(step_node.depth - 1)
                step_node.set_depth(input_node.depth + 1)

                #   4.
                s: Optional[ToolOutput] = None  # the specified cur_step output
                output_dict = input_node.outputs()

                if input_node.node_type == NodeType.OUTPUT:
                    raise Exception(f"Can't connect output {input_node.label} to "
                                    f"input {step_node.label} with tag {input_tag}")

                elif input_node.node_type == NodeType.INPUT:

                    if len(output_dict) == 0:
                        raise Exception(f"Failed when connecting {input_node.label} to {step_node.label}"
                                        f"\n\tImplementation of {input_node.label} incorrectly returns no type for input {input_tag}")
                    if len(output_dict) > 1:
                        raise Exception(f"Implementation of input with {input_node.label} contains "
                                        f"{len(output_dict)} outputs, and we assert exactly 1 output")

                    if len(inp_tag_parts) > 1:
                        Logger.log(f"Step '{step_node.label}' provided too many tags when referencing input: "
                                   f"'{required_step_label}' (superfluous: {', '.join(inp_tag_parts[1:])})",
                                   LogLevel.WARNING)
                        input_tag = required_step_label
                        cur_step.set_input_value(inp_tag, input_tag)

                    # Get the only value in the key
                    s = next(iter(output_dict.values()))

                else:
                    # Now we know we have a STEP node

                    if len(inp_tag_parts) == 1:
                        if len(output_dict) == 1:
                            s = next(iter(output_dict.values()))
                        else:
                            raise Exception(
                                f"The input tag '{input_tag}' requires more information to identify output from '{input_node.label}'")
                    elif inp_tag_parts[1] in output_dict:
                        s = output_dict[inp_tag_parts[1]]
                    elif len(output_dict) == 1:
                        s = next(iter(output_dict.values()))
                        nout = next(iter(output_dict.keys()))
                        input_tag = f"{input_node.label}/{next(iter(output_dict.keys()))}"
                        Logger.log(
                            f"Output '{step_node.label}' did not correctly specify an output of '{input_node.label}', "
                            f"this has been corrected as there was only one output: "
                            f"{inp_tag_parts[1]} → {nout}", LogLevel.WARNING)
                        cur_step.set_input_value(inp_tag, input_tag)

                    else:
                        raise Exception(
                            f"Couldn't uniquely identify input for '{step_node.label}' with output from '{input_node.label}'"
                            f", searching for key: '{inp_tag_parts[1]}' with {', '.join(output_dict.keys())}")

                if s is None:
                    raise Exception("An internal error occurred when determining the connection between "
                                    f"{input_node.label} to {step_node.label} with tag {input_tag}")
                input_type = required_inputs[inp_tag].input_type
                correct_type = input_type.can_receive_from(s.output_type)
                if not correct_type:
                    Logger.log(f"Mismatch of types when joining '{input_tag}' to {step_node.label}.{inp_tag}' "
                               f"({s.output_type.id()} -/→ {required_inputs[inp_tag].input_type.id()})",
                               LogLevel.CRITICAL)
                    Logger.log(f"No action taken to correct type-mismatch of '{input_node.label}' "
                               f"to {step_node.label}/{input_tag}'")

                # # Todo: generate unique identifier of edge (startLabel/tag > endLabel/input)
                # key = f"{input_tag}>{step_node.label}/{inp_tag}"
                col = 'black' if correct_type else 'r'
                work_graph.add_edge(input_node, step_node, type_match=correct_type, color=col)

        for output_node in output_nodes:
            # We did matching when creating outputs, so we can literally just add edges
            req_node = output_node.output.source.split("/")[0]
            if req_node not in labels:
                raise Exception(f"Could not find cur_step {req_node} when stitching output for '{output_node.label}'")
            step_node = labels[req_node]
            output_node.set_depth(step_node.depth + 1)
            work_graph.add_edge(step_node, output_node, color='black')

        end = timer()
        Logger.log("Built DAG ({:2f} s)".format(end - start), LogLevel.INFO)

        return work_graph

    def cwl(self):
        # Let's try to emit CWL
        d = {
            "class": "Workflow",
            "cwlVersion": "v1.0",
            "id": self.label,
            "label": self.label,
            "requirements": [
                {"class": "InlineJavascriptRequirement"}
            ]
        }

        if self.inputs:
            d["inputs"] = {i.label: i.cwl() for i in self.inputs}

        if self.outputs:
            d["outputs"] = {o.label: o.cwl() for o in self.outputs}

        if self.steps:
            d["steps"] = {s.id(): s.cwl() for s in self.steps}

        tools = []
        tools_to_build: Dict[str, Tool] = {s.get_tool().tool(): s.get_tool() for s in self.steps}
        for t in tools_to_build:
            tools.append(tools_to_build[t].cwl())

        inp = {i.id(): i.input_cwl_yml() for i in self.inputs}

        return d, inp, tools

    def wdl(self):

        get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])

        tools = [s.get_tool() for s in self.steps]
        tool_name_to_tool: Dict[str, Tool] = {t.tool().lower(): t for t in tools}
        tool_name_to_alias = {}
        steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in self.steps}

        aliases = set()

        for tool in tool_name_to_tool:
            a = get_alias(tool).upper()
            s = a
            idx = 2
            while s in aliases:
                s = a + idx
                idx += 1
            aliases.add(s)
            tool_name_to_alias[tool] = s

        tab_char = '\t'
        nline_char = '\n'

        imports = '\n'.join([f"import tools/{t}.wdl as {tool_name_to_alias[t.lower()].upper()}" for t in tool_name_to_tool])
        inputs = '\n'.join([f"{tab_char}{i.data_type.primitive()} {i.label}" for i in self.inputs])
        steps = '\n'.join([f"{tab_char}call {tool_name_to_alias[s.get_tool().tool().lower()].upper()}"
                   f".{s.get_tool().tool()} as {s.id()} {{\n{(',' + nline_char).join([2 * tab_char + w for w in s.wdl_map()])}\n}}" for s in self.steps])
        outputs = '\n'.join([f"\t\t{o.data_type.primitive()} {o.label} = {steps_to_alias[o.source.split('/')[0].lower()].lower()}" \
                   f".{o.label}" for o in self.outputs])

        workflow = f"""
{imports}

workflow {self.label} {{
{inputs}
{steps}

    output {{
{outputs}
    }}
}}"""

        tools = [t.wdl() for t in tools]
        inp = {f"{self.label}.{i.id()}": i.data_type.get_value_from_meta(i.meta) for i in self.inputs}

        return workflow, inp, tools


    def _dump_graph(self):
        work_graph = self.__work_graph
        # tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
        tree = json_graph.node_link_data(work_graph, {'link': 'flow', 'source': 'step', 'target': 'target'})
        Logger.log("Workflow Graph: [")
        Logger.log(str(tree))
        Logger.log("] End Workflow Graph")

    def _check_translated(self):
        if self.__work_graph is None:
            raise errors.PipelineTranslatorException(
                'called a translation method before attempting to extract translated components.')

    @staticmethod
    def _translate_step_to_target(step, mapped_inputs):
        return step.translate(mapped_inputs)

    def translate_pipeline_to_json(self):
        json_doc = self._translate_workflow_to_json()
        pretty_json_text = json.dumps(json_doc, indent=2)
        return pretty_json_text

    def _translate_workflow_to_json(self):
        work_graph = self.__work_graph
        Logger.log("Generating JSON description for workflow")
        ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
        type_attr_map = nx.get_edge_attributes(work_graph, 'type')
        workflow = self._describe_workflow(work_graph, self.__input_step, type_attr_map, ctx_attr_map)
        json_doc = {"workflow": workflow}
        Logger.log("Done generating JSON description for workflow")
        return json_doc

    @staticmethod
    def _output_doc_from(step_ctx, step):
        doc = {}
        provides = step_ctx.provides_for(step)
        for output in provides:
            doc[output[Step.STR_ID]] = {
                'type': output[Step.STR_TYPE]
            }

        return doc

    def translate_file(self, pdfile):
        pdfile_path = os.path.abspath(pdfile)
        Logger.log("Using pipeline definition file: " + pdfile_path)

        # Check if the file to translate exists
        if not os.path.isfile(pdfile_path):
            raise ValueError("Specified pipeline definition file does not exist.")

        # All good so lets starts the translation
        with open(pdfile, 'r') as ifile:
            # Validate YAML syntax
            doc = yaml.load(ifile)

        self.__do_translate(doc)

    def translate_string(self, in_string: str):
        """
        Load the YAML string
        :param in_string:
        :return:
        """
        # Converts string into memory stream
        f = io.StringIO(in_string)
        doc = yaml.load(f)
        self.__do_translate(doc)
