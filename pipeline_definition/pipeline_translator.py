import json
import os

import networkx as nx
import yaml
from cerberus import Validator
from networkx.readwrite import json_graph

from pipeline_definition.types.Input_step import InputStep
from pipeline_definition.types.schema import schema
from pipeline_definition.types.step_type import Step
from pipeline_definition.types.type_registry import get_input_factory
from pipeline_definition.types.type_registry import get_step_factory
from pipeline_definition.utils.step_context import StepContext


class PipelineTranslatorException(Exception):
  pass


class PipelineTranslator:
  def __init__(self):
    self.__input_step = None

  @staticmethod
  def __dump_yaml(doc):
    # Diagnostic - what have we got?
    print("YAML DOC [")
    print(yaml.dump(doc))
    print("] END YAML DOC")

  @staticmethod
  def __dump_schema(this_schema):
    print("PDX SCHEMA [")
    # print(schema)
    print(json.dumps(this_schema, indent=4))
    print("] END PDX SCHEMA")

  def validate_schema(self, yaml_doc):
    self.__dump_yaml(yaml_doc)

    sch = schema()
    self.__dump_schema(sch)

    v = Validator(sch)

    validation_success = v.validate(yaml_doc)

    if not validation_success:
      msg = "ERROR! Pipeline definition document validation failed."
      print(msg)
      if v.errors is not None:
        raise ValueError(msg, v.errors)

  @staticmethod
  def build_inputs(inputs):

    input_set = list()

    for input_id, meta in inputs.items():

      if isinstance(meta, str):
        input_type = meta
        meta = dict([(input_type, None)])
      elif isinstance(meta, dict):
        input_type = next(iter(meta.keys()))
      else:
        input_type = None

      print("Processing INPUT: ", input_id, " - ", input_type)

      inp_factory = get_input_factory(input_type)
      if inp_factory is None:
        raise ValueError("No factory registered for input: " + input_type)

      input_obj = inp_factory.build_from(dict([(input_id, meta)]))

      input_set.append(input_obj)

    return input_set

  @staticmethod
  def build_outputs(outputs):
    return None

  @staticmethod
  def build_steps(steps):

    pipeline_steps = list()

    for step in steps:

      step_id = next(iter(step.keys()))
      meta = next(iter(step.values()))

      if isinstance(meta, str):
        step_type = meta
        meta = dict([(step_type, None)])
      elif isinstance(meta, dict):
        step_type = Step.select_type_name_from(meta)
      else:
        step_type = None

      print("Processing STEP: ", step_id, " - ", step_type)

      step_factory = get_step_factory(step_type)
      if step_factory is None:
        raise ValueError("No factory registered for step: " + step_type)

      step_obj = step_factory.build_from(dict([(step_id, meta)]))

      step_obj.validateInputOuputSpec()

      pipeline_steps.append(step_obj)

    return pipeline_steps

  def __create_workflow_graph(self, pipeline_steps, workflow_input_set):
    work_graph = nx.MultiDiGraph()

    # Lets create the input step - the start step of every workflow that produces the workflow inputs as step output
    self.__input_step = InputStep(workflow_input_set)
    pipeline_steps.insert(0, self.__input_step)

    # Now lets put all the steps as node. Then we will establish the edges
    print("Graph Construction: adding nodes for steps.")
    for step in pipeline_steps:
      step_ctx = StepContext(step)
      work_graph.add_node(step, ctx=step_ctx)
      print("Added node for step", step.id(), "[", step.tag(), "]")

    print("Graph Construction: Processing tags.")
    # Now lets put the edges to indicate execution order as indicated by 'tag'
    # Convention is all steps belonging to a tag are exceuted in the seq they have been specified
    # so creating the concept of 'threads' / 'branches'

    # Create a dict of tags vocabulary that have been used in the workflow descriptin
    tag_map = dict()
    for step in pipeline_steps:
      tag = step.tag()
      if tag not in tag_map:
        tag_map[tag] = None

    print("TAG MAP:", str(tag_map))

    # Lets stitch the DAG based on tags assigned to each step

    for step in pipeline_steps:
      print("Processing step ", step.id(), "[", step.tag(), "]")

      # tag specification in step
      stag = step.tag()
      # print("Step specifies tag:", stag)

      if stag == self.__input_step.tag():
        tag_map[stag] = self.__input_step
        continue

      # Have we already seen a step in that thread/tag?
      last_node = tag_map.get(stag)
      tag_map[stag] = step

      if last_node is None:
        last_node = self.__input_step

      work_graph.add_edge(last_node, step, type="branch", tag=stag)
      print("Added edge [", last_node.id(), "]->[", step.id(), "] for tag", stag)

    # Next pass is about input dependency of a step on outputs form other branches and steps
    print("Graph Construction: Processing input dependency")
    for step in pipeline_steps:
      print("Processing step ", step.id(), "[", step.tag(), "]")

      if step == self.__input_step:
        print("Input step never have dependency.")
        continue

      # Lets get the dependency list of the step
      depends = self.dependency_list_of(step)
      if not depends:
        print("Step has no dependency.")
        continue

      for dependency in depends:
        print("Step has dependecny:", dependency)
        self.addDependencyTo(step, dependency, work_graph)

    print("Graph Construction: Polulating context for steps.")
    # Now we have a graph that has edges for flow and dependency
    # Now we need to establish the context of each step
    ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
    type_attr_map = nx.get_edge_attributes(work_graph, 'type')
    self.__populate_step_context(work_graph, self.__input_step, None, ctx_attr_map, type_attr_map)

    return work_graph

  def __populate_step_context(self, work_graph, step, prev_step, ctx_attr_map, type_attr_map):
    print("Populating context for step:", step.id(), "in branch [", step.tag(), "]")
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

    print("Step [", step.id(), "] in branch [", step.tag(), "] has", len(edges), "edges")

    # First pass to process the input dependency
    for edge in edges:
      if edge[0] != step:
        raise RuntimeError("Trouble in edge finding")
      edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
      if edge_type != 'dependency':
        continue
      dependent_on = edge[1]
      dependent_on_ctx = ctx_attr_map[dependent_on]
      print("Step", step.id() + "[" + step.tag() + "] has dependency on step", dependent_on.id(), "[", dependent_on.tag(),
            "]")
      step_ctx.add_dependency_context_from(dependent_on_ctx)

    step_ctx.print()

    # Second pass to process the flow
    for edge in edges:
      if edge[0] != step:
        raise RuntimeError("Trouble in edge finding")
      edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
      if edge_type != 'branch':
        continue
      self.__populate_step_context(work_graph, edge[1], step, ctx_attr_map, type_attr_map)

    return

  def addDependencyTo(self, step, dependency_spec, work_graph):

    if not step:
      return

    if not dependency_spec:
      return

    if not work_graph:
      return

    target_tag = dependency_spec.get('tag')

    # Lets start on the graph by tag
    source_step = self.__input_step
    edge_mapfor_tags = nx.get_edge_attributes(work_graph, 'tag')

    last_step = self.last_step_in_tag(target_tag, source_step, work_graph, edge_mapfor_tags)
    if last_step:
      print(step.id(), "in branch [", step.tag(), "] has input dependency on step", last_step.id(), "in branch [",
            target_tag, "]")
      work_graph.add_edge(step, last_step, type="dependency")

    return

  def last_step_in_tag(self, tag, root, workGraph, edgeMapforTags):
    source_step = root
    edges = nx.edges(workGraph, source_step)
    if not edges:
      # there are no edges so this is the last step
      return source_step

    for edge in edges:
      edge_tag = str(edgeMapforTags[(edge[0], edge[1], 0)])
      if edge_tag == tag:
        return self.last_step_in_tag(tag, edge[1], workGraph, edgeMapforTags)
    return source_step

  @staticmethod
  def dependency_list_of(step):
    if not step:
      return None

    step_requires = step.requires()
    if not step_requires:
      return None

    dependency_list = None
    for requirement in step_requires:
      requirement_name = requirement[Step.STR_ID]
      # print("Process STEP REQUIREMENT:", requirement_name)
      requirement_value = step.providedValueForRequirement(requirement_name)
      # print("Input value:", requirement_value)

      if not requirement_value:
        continue

      dependency_spec = Step.dependency_spec_from(requirement_value)

      if not dependency_list:
        dependency_list = list()
      dependency_list.append(dependency_spec)

    return dependency_list

  @staticmethod
  def __dump_graph(workGraph):
    # tree = json_graph.tree_data(workGraph, self.__root, attrs={'children': 'next', 'id': 'step'})
    tree = json_graph.node_link_data(workGraph, {'link': 'flow', 'source': 'step', 'target': 'target'})
    print("Workflow Graph: [")
    print(tree)
    # jsonDoc = json.dumps(tree, indent=4)
    # print(jsonDoc)

    print("] End Workflow Graph")

  def translate_pipeline(self, pipeline_steps, global_input_set, global_output_set):

    work_graph = self.__create_workflow_graph(pipeline_steps, global_input_set, global_output_set)
    self.__dump_graph(work_graph)

    json_doc = self.translate_workflow_to_json(work_graph)

    pretty_json_text = json.dumps(json_doc, indent=4)

    return pretty_json_text

  def translate_workflow_to_json(self, work_graph):
    print("Generating JSON description for workflow")
    ctx_attr_map = nx.get_node_attributes(work_graph, 'ctx')
    type_attr_map = nx.get_edge_attributes(work_graph, 'type')
    workflow = self.describe_workflow(work_graph, self.__input_step, type_attr_map, ctx_attr_map)
    json_doc = {"workflow": workflow}
    print("Done generating JSON description for workflow")
    return json_doc

  def describe_workflow(self, work_graph, input_step, type_attr_map, ctx_attr_map):

    step_ctx = ctx_attr_map[input_step]
    doc = {'inputs': self.__output_doc_from(step_ctx, input_step)}

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
        self.__populate_description_from(edge[1], branch_desc, work_graph, 1, type_attr_map, ctx_attr_map)
        flow[branch_tag] = {
          'steps': branch_desc
        }

    doc['flow'] = flow

    return doc

  @staticmethod
  def __output_doc_from(step_ctx, step):
    doc = {}
    provides = step_ctx.provides_for(step)
    for output in provides:
      doc[output[Step.STR_ID]] = {
        'type': output[Step.STR_TYPE]
      }

    return doc

  def __populate_description_from(self, step, doc, work_graph, step_order, type_attr_map, ctx_attr_map):

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
    desc['step-outputs'] = self.__output_doc_from(step_ctx, step)

    edges = nx.edges(work_graph, step)
    if not edges:
      return

    for edge in edges:
      if edge[0] != step:
        raise RuntimeError("Trouble in edge finding")
      edge_type = str(type_attr_map[(edge[0], edge[1], 0)])
      if edge_type == 'branch':
        next_step = edge[1]
        self.__populate_description_from(next_step, doc, work_graph, step_order + 1, type_attr_map, ctx_attr_map)

  def translate_yaml_doc(self, yamlDoc):
    # Create a in memory instances for all the inputs, steps and outputs

    inputs = yamlDoc.get("inputs")
    steps = yamlDoc.get("steps")
    outputs = yamlDoc.get("outputs")

    if inputs is None:
      raise ValueError("No input?")

    # if outputs is None:
    #    raise ValueError("No output?")

    workflow_input_set = self.build_inputs(inputs)
    workflow_output_set = self.build_outputs(outputs)
    pipeline_steps = self.build_steps(steps)

    # Now translate the workflow steps
    doc = self.translate_pipeline(pipeline_steps, workflow_input_set, workflow_output_set)

    return doc

  def translate(self, pdfile, outfile=None, overwrite_outfile=False):
    pdfile_path = os.path.abspath(pdfile)
    print("Using PD file: " + pdfile_path)

    # Check if the file to translate exists
    if not os.path.isfile(pdfile_path):
      raise ValueError("Specified PD file does not exist.")

    # The output file name, if not explicitly specified, a convention is used
    if outfile is None:
      outfile = pdfile + ".pdx"
    outfile_path = os.path.abspath(outfile)
    print("Using Output file: " + outfile_path)

    # If the output file path already exists as directory?
    if os.path.isdir(outfile_path):
      raise ValueError("Directory already exists with the name of outfile.")

    # Otherwise can we overwrite?
    if overwrite_outfile is False:
      if os.path.isfile(outfile_path):
        raise ValueError("Outfile already exists.")

    # All good so lets starts the translation
    with open(pdfile, 'r') as ifile:
      # Validate YAML syntax
      doc = yaml.load(ifile)

      # Do schema validation
      self.validate_schema(doc)

      # Doc is OK, lets translate
      output_doc = self.translate_yaml_doc(doc)
      print("Translation Output: [")
      print(output_doc)
      print("]")

      # Save it
      # with open( outfile_path, "w" ) as ofile:
      #    ofile.write( tDoc )
