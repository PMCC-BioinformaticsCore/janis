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
from typing import List, Dict, Any, Optional, Type
from timeit import default_timer as timer


import networkx as nx
from cerberus import Validator
from networkx.readwrite import json_graph

import pipeline_definition.types.schema as wehi_schema
import pipeline_definition.utils.errors as errors
from pipeline_definition.types.common_data_types import String, Number, Boolean, Array
from pipeline_definition.types.data_types import DataType
from pipeline_definition.types.output import Output, OutputNode

from pipeline_definition.utils.logger import Logger, LogLevel

from pipeline_definition.graph.node import Node, NodeType
from pipeline_definition.types.input import InputNode
from pipeline_definition.utils.yaml_utils import str_presenter
from pipeline_definition.types.type_registry import get_tool, get_type
# from pipeline_definition.types.type_registry import get_input_factory
from pipeline_definition.types.step import Step, ToolOutput, StepNode
from pipeline_definition.types.input import Input  #, InputFactory, InputType

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
    def __init__(self):
        Logger.log("Creating PipelineTranslator", LogLevel.CRITICAL)

        self.__input_step: InputNode = None
        self.__work_graph: nx.Graph = None
        self.__labels_map: Dict[str, Node] = {}

        self.__inputs: List[Input] = []
        self.__steps: List[Step] = []
        self.__outputs: List[Output] = []

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

        workflow_input_set: List[Input] = self._build_inputs(inputs)
        pipeline_steps: List[Step] = self._build_steps(steps)
        workflow_output_set: List[Output] = self._build_outputs(outputs, workflow_input_set, pipeline_steps)

        self.__work_graph = self._create_workflow_graph_and_outputs(pipeline_steps, workflow_input_set, workflow_output_set)
        self.draw_graph()
        self._dump_graph()

        # # Now translate the workflow steps

    def draw_graph(self):
        import matplotlib.pyplot as plt

        default_color = 'blacl'

        G = self.__work_graph
        edges_attributes = [G.edges[e] for e in G.edges]
        edge_colors = [x["color"] if 'color' in x else default_color for x in edges_attributes]
        node_colors = [NodeType.to_col(x.node_type) for x in G.nodes]

        nx.draw(G, edge_color=edge_colors, node_color=node_colors, with_labels=True)
        plt.show()

    def _build_inputs(self, inputs: Dict[str, Any]) -> List[Input]:
        """
        Build a list of inputs from the inputs dictionary
        :param inputs: Dict[str, Any], where Any = str | { $Type: InputData }
        :return: List[Input] >> All the inputs in the workflow
        """
        input_set: List[Input] = []
        Logger.log("Building inputs", LogLevel.DEBUG)

        for input_id, meta in inputs.items():           # StepLabel, StepProperties/Data

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
                raise Exception(f"The step '{req_node}' does not contain any outputs to connect to the output: '{output_id}'")

            if len(provides) == 1:
                so = provides[next(iter(provides))]
                if (len(inp_tag_parts)) == 1:
                    new_output_source = f"{output_source}/{so.tag}"
                    Logger.log(f"Output '{output_id}' should fully specify step/outputname: "
                               f"{output_source} → {new_output_source}", LogLevel.WARNING)
                    output_source = new_output_source
                elif inp_tag_parts[1] not in provides:
                    new_output_source = f"{output_source}/{so.tag}"
                    Logger.log(f"Output '{output_id}' did not correctly specify an output of '{req_node}', "
                               f"this has been corrected as the step '{req_node}' only contained one ouput: "
                               f"{output_source} → {new_output_source}", LogLevel.WARNING)

                output = Output(output_id, output_source, so.output_type)
                outputs_ar.append(output)

            else:
                if len(inp_tag_parts) == 1:
                    raise Exception(f"Output '{output_id}' ({output_source}) did not fully specify step/outputname")
                if inp_tag_parts[1] not in provides:
                    raise Exception(f"The output '{output_id}' could not find an output called '{inp_tag_parts[1]}' in '{req_node}'")
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
        # TODO: Remove type inferencing, and devise strong referencing for step inputs.
        #       For this, we'll need to define a specific format (probably close to CWL)
        #       like: - { tool: type/version, ...inputs }

        pipeline_steps: List[Step] = []

        for step_id, meta in steps.items():

            step_type: Optional[str]

            if isinstance(meta, str):       # I don't think we'll allow this to happen, because it still
                step_type = str(meta)       # implies type inferencing to some point, or maybe all outputs go to this
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
                raise ValueError("No factory registered for tool: " + step_type)

            # step_obj = step_factory.build_from(step_id, meta)
            step_obj = Step(step_id, tool_type(), meta)
            pipeline_steps.append(step_obj)

        return pipeline_steps

    def _create_workflow_graph_and_outputs(self, pipeline_steps: List[Step], workflow_inputs: List[Input], workflow_outputs: List[Output]) -> nx.Graph:
        """
        Take the inputs, {outputs} and steps and stitch it all together
        :param pipeline_steps: List[Step]
        :param workflow_inputs: List[Input]
        :return: a networx graph
        """
        start = timer()
        Logger.log("Building DAG...", LogLevel.INFO)
        work_graph = nx.MultiDiGraph()                  # Use the DAG from networx

        # Lets create the input step - the start step of every workflow that produces the workflow inputs as step output
        # self.__input_step = InputNode(workflow_inputs)
        # pipeline_steps.insert(0, self.__input_step)

        input_nodes = [InputNode(inp) for inp in workflow_inputs]
        work_graph.add_nodes_from(input_nodes)

        step_nodes = [StepNode(step) for step in pipeline_steps]
        work_graph.add_nodes_from(step_nodes)

        output_nodes = [OutputNode(output) for output in workflow_outputs]
        work_graph.add_nodes_from(output_nodes)

        labels = { n.label: n for n in work_graph.nodes }
        Logger.log(f"Node labels: {labels}")

        # Todo: Move the matching and type checking logic to the Input creation

        for step_node in step_nodes:
            # Create edges
            required_inputs = step_node.step.requires()
            for inp_tag in required_inputs:
                inp = step_node.step.input_value(inp_tag)
                # Get the correct label, to build the acylic graph, we really only need to first section when split by /
                inp_tag_parts = inp.split("/")
                required_step_label = inp_tag_parts[0]

                input_node: Node = labels[required_step_label]

                # Todo: generate unique identifier of edge (startLabel/tag > endLabel/input) + check if types match
                key = f"{inp}>{step_node.label}/{inp_tag}"

                # Check types
                s: Optional[ToolOutput] = None    # the specified step output
                if input_node.node_type == NodeType.OUTPUT:
                    raise Exception(f"Can't connect output {input_node.label} to "
                                    f"input {step_node.label} with tag {inp}")
                elif input_node.node_type == NodeType.INPUT:
                    output_dict = input_node.outputs()
                    # Get the only value in the key
                    if len(output_dict) == 0:
                        raise Exception(f"Failed when connecting {input_node.label} to {step_node.label}"
                            f"\n\tImplementation of {input_node.label} incorrectly returns no type for input {inp}")
                    if len(output_dict) > 1:
                        raise Exception(f"Implementation of input with {input_node.label} contains "
                                        f"{len(output_dict)} outputs, and we assert exactly 1 output")

                    s = next(iter(output_dict.values()))

                else:
                    # Now we know we have a STEP node
                    output_dict = input_node.outputs()

                    if len(inp_tag_parts) == 1:
                        if len(output_dict) == 1:
                            s = next(iter(output_dict.values()))
                        else:
                            raise Exception(f"The input tag '{inp}' requires more information to identify output from '{input_node.label}'")
                    elif inp_tag_parts[1] in output_dict:
                        s = output_dict[inp_tag_parts[1]]
                    else:
                        raise Exception(f"Couldn't uniquely identify output from '{input_node.label}' for '{step_node.label}'"
                                        f", searching for key: '{inp_tag_parts[1]}'")

                if s is None:
                    raise Exception("An internal error occurred when determining the connection between "
                                    f"{input_node.label} to {step_node.label} with tag {inp}")

                correct_type = s.output_type.can_receive_from(required_inputs[inp_tag].input_type)
                if not correct_type:
                    Logger.log(f"Mismatch of types between {input_node.label} ({s.output_type}) to "
                               f"{step_node.label} with tag {inp} ({required_inputs[inp_tag].input_type})",
                               LogLevel.CRITICAL)
                col = 'black' if correct_type else 'r'
                work_graph.add_edge(input_node, step_node, type_match=correct_type, color=col)

        for output_node in output_nodes:
            # We did matching when creating outputs, so we can literally just add edges
            req_node = output_node.output.source.split("/")[0]
            if req_node not in labels:
                raise Exception(f"Could not find step {req_node} when stitching output for '{output_node.label}'")
            step_node = labels[req_node]
            work_graph.add_edge(step_node, output_node, color='black')

        end = timer()
        Logger.log("Built DAG ({:2f} s)".format(end - start), LogLevel.INFO)

        return work_graph

    def _dump_graph(self):
        work_graph = self.__work_graph
        # tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
        tree = json_graph.node_link_data(work_graph, {'link': 'flow', 'source': 'step', 'target': 'target'})
        Logger.log("Workflow Graph: [")
        Logger.log(str(tree))
        Logger.log("] End Workflow Graph")

    def _check_translated(self):
        if self.__work_graph is None:
            raise errors.PipelineTranslatorException('called a translation method before attempting to extract translated components.')

    def input(self, resolve=False):
        self._check_translated()

        input_items = self.__input_step.inputs()

        s = dict()
        for inp in input_items:
            if resolve:
                inp.resolve()
            s.update(inp.translate_for_input())

        inp = {'inputs': s}

        return yaml.dump(inp, default_flow_style=False)

    def pipeline(self):
        self._check_translated()

        work_graph = self.__work_graph
        input_step = self.__input_step

        ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
        type_attr_map = nx.get_edge_attributes(work_graph, 'type')
        return self._start_workflow_translation(work_graph, input_step, type_attr_map, ctx_attr_map)

    def _start_workflow_translation(self, work_graph, input_step, type_attr_map, ctx_attr_map):

        preamble = yaml.load("""
class: Workflow
cwlVersion: v1.0

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var rename_trim_file = function() {
      if ( self == null ) {
        return null;
      } else {
        var xx = self.basename.split('.');
        var id = xx.indexOf('fastq');
        xx.splice(id, 1);
        return xx.join('.');
      }
    };
- class: ScatterFeatureRequirement
- class: StepInputExpressionRequirement
- class: SchemaDefRequirement
  types:
  - $import: ../tools/src/tools/trimmomatic-end_mode.yml
  - $import: ../tools/src/tools/trimmomatic-sliding_window.yml
  - $import: ../tools/src/tools/trimmomatic-phred.yml
  - $import: ../tools/src/tools/trimmomatic-illumina_clipping.yml
  - $import: ../tools/src/tools/trimmomatic-max_info.yml
""")

        inputs = self._translate_step_to_target(input_step, None)

        edges = nx.edges(work_graph, input_step)
        if not edges:
            return

        steps_xlate = dict()

        for edge in edges:
            if edge[0] != input_step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type == 'branch':
                self._recurse_graph_and_translate(edge[1], work_graph, 1, type_attr_map, ctx_attr_map, steps_xlate)

        outputs = self._find_step_outputs(steps_xlate)

        xlate = dict()
        xlate.update(preamble)
        xlate.update(inputs)
        xlate.update({'steps': steps_xlate})
        xlate.update(outputs)

        return yaml.dump(xlate, default_flow_style=False)

    @staticmethod
    def _find_step_outputs(steps):
        # Find the outputs from a dictionary of CWL step and generate the CWL output entry
        outputs = dict()
        for (step_id, step) in steps.items():
            try:
                outs = step['out']
                for out in outs:
                    outputs.update({f'{step_id}_{out}': {'type': 'File', 'outputSource': f'{step_id}/{out}'}})
            except KeyError:
                print(f'Step: {step_id} has no outputs.')

        return {'outputs': outputs}

    def _recurse_graph_and_translate(self, step, work_graph, step_order, type_attr_map, ctx_attr_map, steps_xlate):

        mapped_inputs = []
        step_inputs = step.requires()

        for step_input in step_inputs:
            step_ctx = ctx_attr_map[step]
            mapping = step_ctx.map_input_for_translation(step_input)

            mi = MappedInput(step_inputs, mapping, step_input)
            mapped_inputs.append(mi)

        s = self._translate_step_to_target(step, mapped_inputs)
        steps_xlate.update(s)

        edges = nx.edges(work_graph, step)
        if not edges:
            return

        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type == 'branch':
                next_step = edge[1]
                self._recurse_graph_and_translate(next_step, work_graph, step_order + 1, type_attr_map, ctx_attr_map,
                                                  steps_xlate)

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

    def _describe_workflow(self, work_graph, input_step, type_attr_map, ctx_attr_map):

        step_ctx = ctx_attr_map[input_step]
        doc = {'inputs': self._output_doc_from(step_ctx, input_step)}

        flow = {}
        edges = nx.edges(work_graph, input_step)
        if not edges:
            return

        for edge in edges:
            if edge[0] != input_step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type == 'branch':
                branch_tag = edge[1].tag()
                branch_desc = {}
                self._populate_description_from(edge[1], branch_desc, work_graph, 1, type_attr_map, ctx_attr_map)
                flow[branch_tag] = {
                    'steps': branch_desc
                }

        doc['flow'] = flow

        return doc

    def _populate_description_from(self, step, doc, work_graph, step_order, type_attr_map, ctx_attr_map):

        desc = {
            'step': step.id(),
            'type': step.type()
        }

        # desc['order'] = stepOrder
        doc[step_order] = desc

        step_inputs = step.requires()
        if step_inputs:
            inputs_doc = {}

            for step_input in step_inputs:
                input_id = step_input[Step.STR_ID]
                input_type = step_input[Step.STR_TYPE]
                idoc = {
                    'type': input_type
                }

                step_ctx = ctx_attr_map[step]
                mapping = step_ctx.map_input(step_input)
                idoc['mapping'] = mapping

                inputs_doc[input_id] = idoc

            desc['step-inputs'] = inputs_doc

        step_ctx = ctx_attr_map[step]
        desc['step-outputs'] = self._output_doc_from(step_ctx, step)

        edges = nx.edges(work_graph, step)
        if not edges:
            return

        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type == 'branch':
                next_step = edge[1]
                self._populate_description_from(next_step, doc, work_graph, step_order + 1, type_attr_map, ctx_attr_map)

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
