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

import json
import os
import io
from typing import List, Dict, Any, Tuple

import networkx as nx
import yaml
from cerberus import Validator
from networkx.readwrite import json_graph

from pipeline_definition.types.input_step import InputStep
from pipeline_definition.types.input_type import InputFactory, Input, InputType
import pipeline_definition.types.schema as WEHI_Schema
from pipeline_definition.types.step_type import Step, DependencySpec
from pipeline_definition.types.type_registry import get_input_factory
from pipeline_definition.types.type_registry import get_step_factory
from pipeline_definition.utils.errors import PipelineTranslatorException
from pipeline_definition.utils.step_context import StepContext
from pipeline_definition.utils.yaml_utils import str_presenter

yaml.add_representer(str, str_presenter)


class MappedInput:
    def __init__(self, inputs: List[InputType], candidates, step_input: InputType):
        self.inputs = inputs
        self.candidates = candidates
        self.input_type = step_input.type_name()

    def __repr__(self):
        return f'inputs={self.inputs}: candidates={self.candidates}'


class PipelineTranslator:
    def __init__(self, debug: bool = False):
        self.__input_step: InputStep = None
        self.__debug: bool = debug
        self.__work_graph: nx.Graph = None

    def _debug_print(self, *args) -> None:
        """
        It would be great to potentially move this to a singleton and
        having levels of errors that this class could auto filter out.
        Potentially something like:
            - None
            - Critical (RED)
            - Information (WHITE)
            - Debug (GRAY)
        :param args:
        :return: None
        """
        if self.__debug:
            print(args)

    def _dump_as_yaml(self, doc: Dict[str, Any], indent: int = 2, prefix="YAML DOC") -> None:
        """
        Dump the provided dictionary as a YAML document
        :param doc: Dict[str, Any]
        """
        # Diagnostic - what have we got?
        self._debug_print(f"{prefix} [")
        self._debug_print(yaml.dump(doc, indent=indent))
        self._debug_print(f"] END {prefix}")

    def __do_translate(self, doc: Dict[str, Any]) -> None:
        """
        Convert dictionary (of WEHI pipeline definition) into DAG
        :param doc: Dict[str, str]
        :return: N/A
        """
        # Do schema validation
        self.validate_schema(doc)

        # Now convert the YAML doc?
        self._translate_yaml_doc(doc)
        self._debug_print("Translation Output: [")
        self._debug_print(self.__work_graph)
        self._debug_print("]")

    def validate_schema(self, yaml_doc: Dict[str, Any]) -> bool:
        """
        Validate the yaml_doc against the WEHI spec,
        throws exception if cerberus gives errors.
        :param yaml_doc:
        :return:
        """
        self._dump_as_yaml(yaml_doc)

        sch = WEHI_Schema.workflow_schema()
        self._dump_as_yaml(sch, indent=2, prefix="PDX SCHEMA")

        v: Validator = Validator(sch)

        validation_success: bool = v.validate(yaml_doc)

        if validation_success:
            self._debug_print("Schema is Valid")
        else:
            msg = "ERROR! Pipeline definition document validation failed."
            self._debug_print(msg)
            if v.errors is not None:
                raise ValueError(msg, v.errors)
        return validation_success

    def _translate_yaml_doc(self, doc: Dict[str, Any]) -> None:
        """
        Pulls out the inputs, outputs, steps and calls to create workflow_graph
        :param doc: Dict[str, Any]
        """
        # Create a in memory instances for all the inputs, steps and outputs

        inputs: Dict[str, Any] = doc.get(WEHI_Schema.KEYS.INPUTS)
        steps: List[Dict[str, Any]] = doc.get(WEHI_Schema.KEYS.STEPS)
        outputs: Dict[str, Any] = doc.get(WEHI_Schema.KEYS.OUTPUTS)

        if inputs is None:
            # This is probably okay, depends on whether the the steps reference anything
            # Probably better to be a warning which we can log at the end
            raise ValueError("No input?")

        if steps is None or not steps:
            raise ValueError("There are no steps in the diagram")

        # if outputs is None:
        #    raise ValueError("No output?")

        workflow_input_set: List[Input] = self._build_inputs(inputs)
        # workflow_output_set: List[Output] = self._build_outputs(outputs)
        pipeline_steps: List[Step] = self._build_steps(steps)

        self.__work_graph = self._create_workflow_graph(pipeline_steps, workflow_input_set)
        self._dump_graph()

        # # Now translate the workflow steps

    def _build_inputs(self, inputs: Dict[str, Any]) -> List[Input]:
        """
        Build a list of inputs from the inputs dictionary
        :param inputs: Dict[str, Any], where Any = str | { $Type: InputData }
        :return: List[Input] >> All the inputs in the workflow
        """
        input_set: List[Input] = []

        for input_id, meta in inputs.items():           # StepLabel, StepProperties/Data

            input_type: str = None                      # Declare this here to improve clarity inside loop
            if isinstance(meta, str):
                # Our step is just a string, ie, it's the tool
                input_type = meta
                meta = dict([(input_type, None)])

            elif isinstance(meta, dict):
                # mfranklin: I don't think this is extremely clear, I'd probably propose a "type" field, ie:
                # | input_type = meta["type"]
                input_type = next(iter(meta.keys()))
            else:
                raise ValueError(f"Rerecognised type ({type(meta)}) for input")

            self._debug_print(f"Processing INPUT: {input_id} - {input_type}")

            inp_factory: InputFactory = get_input_factory(input_type)
            if inp_factory is None:
                raise ValueError("No factory registered for input: " + input_type)

            # Build an Input from this data
            input_obj = inp_factory.build_from({input_id: meta}, self.__debug)
            input_set.append(input_obj)

        return input_set

    @staticmethod
    def _build_outputs(outputs):
        return None

    def _build_steps(self, steps: List[Dict]) -> List[Step]:
        """
        Build list of steps (unconnected) purely from the workflow
        :param steps:
        :return:
        """
        # TODO: Remove type inferencing, and devise strong referencing for step inputs.
        #       For this, we'll need to define a specific format (probably close to CWL)
        #       like: - { tool: type/version, ...inputs }

        pipeline_steps: List[Step] = []

        for step in steps:

            step_id: str = next(iter(step.keys()))
            meta = next(iter(step.values()))            # First key in the dictionary

            if isinstance(meta, str):
                # With type inferencing, we'll remove this step
                step_type = meta
                meta = dict([(step_type, None)])
            elif isinstance(meta, dict):
                step_type = Step.select_type_name_from(meta)
            else:
                step_type = None

            self._debug_print(f"Processing STEP: {step_id} - {step_type}")

            step_factory = get_step_factory(step_type)
            if step_factory is None:
                raise ValueError("No factory registered for step: " + step_type)

            step_obj = step_factory.build_from({step_id: meta}, debug=self.__debug)
            pipeline_steps.append(step_obj)

        return pipeline_steps

    def _create_workflow_graph(self, pipeline_steps: List[Step], workflow_input_set: List[Input]) -> nx.Graph:
        """
        Take the inputs, {outputs} and steps and stitch it all together
        :param pipeline_steps: List[Step]
        :param workflow_input_set: List[Input]
        :return: a networx graph
        """
        work_graph = nx.MultiDiGraph()                  # Use the DAG from networx

        # Lets create the input step - the start step of every workflow that produces the workflow inputs as step output
        self.__input_step = InputStep(workflow_input_set)
        pipeline_steps.insert(0, self.__input_step)

        # Now lets put all the steps as node. Then we will establish the edges
        self._debug_print("Graph Construction: adding nodes for steps.")
        for step in pipeline_steps:
            step_ctx = StepContext(step)
            work_graph.add_node(step, ctx=step_ctx)
            self._debug_print("Added node for step", step.id(), "[", step.tag(), "]")

        self._debug_print("Graph Construction: Processing tags.")

        # Now let's put the edges to indicate execution order as indicated by 'tag'
        # Convention is all steps belonging to a tag are exceuted in the seq they have been specified
        # so creating the concept of 'threads' / 'branches.
        # Create a dict of tags vocabulary that have been used in the workflow description
        tag_map: Dict[str, Any] = { step.tag(): None for step in pipeline_steps }
        # for step in pipeline_steps:
        #     tag = step.tag()
        #     if tag not in tag_map:
        #         tag_map[tag] = None

        self._debug_print("TAG MAP:", str(tag_map))

        # Lets stitch the DAG based on tags assigned to each step

        for step in pipeline_steps:
            self._debug_print("Processing step ", step.id(), "[", step.tag(), "]")

            # tag specification in step
            stag: str = step.tag()
            # self.__debug_print("Step specifies tag:", stag)

            if stag == self.__input_step.tag():
                tag_map[stag] = self.__input_step
                continue

            # Have we already seen a step in that thread/tag?
            last_node = tag_map.get(stag)
            tag_map[stag] = step

            if last_node is None:
                last_node = self.__input_step

            work_graph.add_edge(last_node, step, type="branch", tag=stag)
            self._debug_print("Added edge [", last_node.id(), "]->[", step.id(), "] for tag", stag)

        # Next pass is about input dependency of a step on outputs form other branches and steps
        self._debug_print("Graph Construction: Processing input dependency")
        for step in pipeline_steps:
            self._debug_print("Processing step ", step.id(), "[", step.tag(), "]")

            if step.tag() == self.__input_step.tag():
                self._debug_print("Input step never have dependency.")
                continue

            # Lets get the dependency list of the step
            depends = self.get_dependecy_spec(step)
            if not depends:
                self._debug_print("Step has no dependency.")
                continue

            for dependency in depends:
                self._debug_print("Step has dependecny:", dependency)
                self._add_dependency_to(step, dependency, work_graph)

        self._debug_print("Graph Construction: Polulating context for steps.")
        # Now we have a graph that has edges for flow and dependency
        # Now we need to establish the context of each step
        ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
        type_attr_map = nx.get_edge_attributes(work_graph, 'type')
        self._populate_step_context(work_graph, self.__input_step, None, ctx_attr_map, type_attr_map)

        return work_graph

    def _populate_step_context(self, work_graph, step, prev_step, ctx_attr_map, type_attr_map):
        self._debug_print("Populating context for step:", step.id(), "in branch [", step.tag(), "]")
        step_ctx = ctx_attr_map.get(step)
        if not step_ctx:
            raise RuntimeError("Missing step context in graph. Graph integrity fail.")

        if prev_step:
            prev_ctx = ctx_attr_map.get(prev_step)
            if not prev_ctx:
                raise RuntimeError("Missing step context in graph. Graph integrity fail.")

            if prev_step.tag() != InputStep.input_step_tag_name() and step.tag() != prev_step.tag():
                raise RuntimeError("Branch tag mismatch during context population.")

            step_ctx.inherit_context_of_branch(prev_ctx)

        edges = nx.edges(work_graph, step)
        if not edges:
            return

        self._debug_print("Step [", step.id(), "] in branch [", step.tag(), "] has", len(edges), "edges")

        # First pass to process the input dependency
        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type != 'dependency':
                continue
            dependent_on = edge[1]
            dependent_on_ctx = ctx_attr_map[dependent_on]
            self._debug_print("Step", step.id() + "[" + step.tag() + "] has dependency on step", dependent_on.id(), "[",
                              dependent_on.tag(),
                              "]")
            step_ctx.add_dependency_context_from(dependent_on_ctx)

        if self.__debug:
            step_ctx.print()

        # Second pass to process the flow
        for edge in edges:
            if edge[0] != step:
                raise RuntimeError("Trouble in edge finding")
            edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
            if edge_type != 'branch':
                continue
            self._populate_step_context(work_graph, edge[1], step, ctx_attr_map, type_attr_map)

        return

    def _add_dependency_to(self, step: Step, dependency_spec: DependencySpec, work_graph: nx.MultiDiGraph) -> None:
        """
        Work out where step connects to
        :param step: Step
        :param dependency_spec: DependencySpec
        :param work_graph: nx.MultiDiGraph
        :return: None
        """

        if not step:
            raise Exception("No step provided to dependency")

        if not dependency_spec:
            raise Exception("No dependency spec provided for add_dependency")

        if not work_graph:
            raise Exception("No work_graph provided for add_dependency")

        target_tag: str = dependency_spec.tag

        # Lets start on the graph by tag
        source_step: InputStep = self.__input_step
        # dictionary for multidigraph is keyed by (u, v, key)
        # https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.classes.function.get_edge_attributes.html
        # TODO: We can create a type that represents the (u, v, key) data structure
        edge_mapfor_tags: Dict[Tuple, str] = nx.get_edge_attributes(work_graph, 'tag')

        last_step: Step = self._last_step_in_tag(target_tag, source_step, work_graph, edge_mapfor_tags)
        if last_step:
            self._debug_print(f"({step.id()}) in branch [{step.tag()}] has input " # automatically joins lines
                              f"dependency on step ({last_step.id()}) in branch [{target_tag}]")
            work_graph.add_edge(step, last_step, type="dependency")

        return

    def _last_step_in_tag(self, tag: str, source_step: InputStep, work_graph: nx.MultiDiGraph, edge_mapfor_tags: Dict[Tuple, str]) -> Step:
        """
        Recursive function to determine the latest step in the branch, specified by $tag.
        :param tag: Branch tag (eg: #cancer)
        :param source_step: InputStep
        :param work_graph: nx.MultiDiGraph
        :param edge_mapfor_tags: Dict[(u, v, key), str]
        :return: Step
        """
        edges: List[Tuple] = nx.edges(work_graph, source_step)
        if not edges:
            # there are no edges so this is the last step
            return source_step

        for edge in edges:
            # dictionary for multidigraph is keyed by (u, v, key)
            # https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.classes.function.get_edge_attributes.html
            edge_tag = str(edge_mapfor_tags[(edge[0], edge[1], 0)])
            if edge_tag == tag:
                return self._last_step_in_tag(tag, edge[1], work_graph, edge_mapfor_tags)
        return source_step

    @staticmethod
    def get_dependecy_spec(step: Step) -> List[DependencySpec]:
        """
        Get a list of dependencies based on the tag from the current step
        :param step: Step
        :return: List[DependencySpec]
        """
        if not step:
            return []

        step_requires: List[InputType] = step.requires()
        if not step_requires:
            return []

        if not isinstance(step_requires, list):
            raise Exception(f"Implementation of .requires() for {type(step)} "
                            f"is not returning a List[InputType] but a ({type(step_requires)})")

        dependency_list: List[DependencySpec] = []

        for requirement in step_requires:
            requirement_name = requirement.type_name()
            requirement_value = step.provided_value_for_requirement(requirement)

            if not requirement_value:
                continue

            dependency_spec = Step.dependency_spec_from(requirement_value)
            dependency_list.append(dependency_spec)

        return dependency_list

    def _dump_graph(self):
        work_graph = self.__work_graph
        # tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
        tree = json_graph.node_link_data(work_graph, {'link': 'flow', 'source': 'step', 'target': 'target'})
        self._debug_print("Workflow Graph: [")
        self._debug_print(tree)
        # jsonDoc = json.dumps(tree, indent=2)
        # self.__debug_print(jsonDoc)

        self._debug_print("] End Workflow Graph")

    def _check_translated(self):
        if self.__work_graph is None:
            raise PipelineTranslatorException(
                'call a translation method before attempting to extract translated components.')

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
        self._debug_print("Generating JSON description for workflow")
        ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
        type_attr_map = nx.get_edge_attributes(work_graph, 'type')
        workflow = self._describe_workflow(work_graph, self.__input_step, type_attr_map, ctx_attr_map)
        json_doc = {"workflow": workflow}
        self._debug_print("Done generating JSON description for workflow")
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
        self._debug_print("Using pipeline definition file: " + pdfile_path)

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
