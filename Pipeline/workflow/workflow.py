import networkx as nx
from typing import Dict, List, Tuple, Set, Optional, Any

from Pipeline.utils import first_value
from Pipeline.utils.logger import Logger, LogLevel

from Pipeline.graph.stepinput import StepInput
from Pipeline.graph.node import Node, NodeTypes, NodeAndTag, layout_nodes2

from Pipeline.types.common_data_types import Array
from Pipeline.types.data_types import DataType, NativeTypes
from Pipeline.utils.descriptor import BaseDescriptor, GetOnlyDescriptor

from Pipeline.workflow.step import Step, StepNode
from Pipeline.workflow.input import Input, InputNode
from Pipeline.workflow.output import Output, OutputNode
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput, ToolTypes, ToolType
from Pipeline.utils.errors import DuplicateLabelIdentifier, InvalidNodeIdentifier, NodeNotFound, InvalidStepsException, \
    InvalidInputsException

import cwlgen.cwlgen as cwl


class Workflow(Tool):

    identifier = BaseDescriptor()
    label = BaseDescriptor()
    doc = BaseDescriptor()

    def __init__(self, identifier: str, label: str = None, doc: Optional[str] = None):
        """
        Initialise the worklow
        :param identifier: uniquely identifies the workflow
        :param label: a label that the engine may use to represent the workflow (should be unique)
        :param doc: Documentation of the workflow
        """
        Logger.log(f"Creating workflow with identifier: '{identifier}'")

        self.identifier = identifier
        self.label = label
        self.doc = doc

        self._nodes: Dict[str, Node] = {}       # Look up a node by its identifier

        self._inputs: List[InputNode] = []      # InputNodes
        self._steps: List[StepNode] = []        # StepNodes
        self._outputs: List[OutputNode] = []    # OutputNodes

        self.graph: nx.MultiDiGraph = nx.MultiDiGraph()     # Realistically this isn't really used except for an image

        # Flags for different requirements that a workflow might need
        self.has_scatter = False
        self.has_subworkflow = False
        self.has_multiple_inputs = False

    def id(self) -> str:
        """
        Returns the identifier of the workflow
        """
        return self.identifier

    @classmethod
    def type(cls) -> ToolType:
        """
        Returns the ToolType (Workflow | CommandTool | ExpressionTool)
        """
        return ToolTypes.Workflow

    def inputs(self) -> List[ToolInput]:
        """
        List of ToolInputs of the workflow, we can toss out most of the metadata
        about positioning, prefixes, etc that the ToolInput class uses
        """
        return [ToolInput(i.id(), i.input.data_type) for i in self._inputs]

    def outputs(self) -> List[ToolOutput]:
        """
        Similar to inputs, return a list of ToolOutputs of the workflow
        """
        return [ToolOutput(o.id(), o.output.data_type, o.output.doc)
                for o in self._outputs if o.output.data_type is not None]

    def draw_graph(self):
        """
        Use matplotlib to draw the a diagram of the workflow (currently neglects subworkflows).
        This will pause the execution of the program until the graph is closed.
        """
        import matplotlib.pyplot as plt

        default_color = 'black'

        G = self.graph
        edges_attributes = [G.edges[e] for e in G.edges]
        edge_colors = [x["color"] if 'color' in x else default_color for x in edges_attributes]
        node_colors = [NodeTypes.to_col(x.node_type) for x in G.nodes]

        pos = layout_nodes2(list(G.nodes))  # Manual layout engine of the graph, need to make this better

        for n in pos:
            G.node[n]['pos'] = pos[n]

        nx.draw(G, pos=pos, edge_color=edge_colors, node_color=node_colors, with_labels=True)
        plt.show()

    def _add_items(self, *args) -> List[Any]:
        """
        Supports adding (Input | Step | Output | Tuple = construct edge)
        """
        mixed: List[Any] = args

        if isinstance(args[0], list):
            if len(args) > 1:
                raise Exception("Invalid type of arguments")
            mixed = args[0]

        return [self._add_item(i) for i in mixed]

    def _add_item(self, item):
        """
        Generic add whatever into the workflow: (Input | Step | Output | Tuple = construct edge)
        :return:
        """
        if isinstance(item, Input):
            return self._add_inputs([item])[0]
        elif isinstance(item, Step):
            return self._add_steps([item])[0]
        elif isinstance(item, Output):
            return self._add_outputs([item])[0]
        elif isinstance(item, tuple):
            return self.add_edge(item[0], item[1])
        else:
            raise Exception(f"Unexpected type '{type(item)}' passed to 'add_nodes'")

    def _add_inputs(self, inputs: List[Input]) -> List[InputNode]:
        """
        Add a list of inputs into the graph, does nothing except add it into the graph
        """
        ins = []
        for inp in inputs:
            Logger.log(f"Adding input '{inp.id()}' to '{self.identifier}'")
            node: InputNode = InputNode(inp)
            ins.append(node)
            self._add_node(node)
            self._inputs.append(node)
        return ins

    def _add_outputs(self, outputs: List[Output]) -> List[OutputNode]:
        """
        Adds a list of output to the graph, and returns the created OutputNodes
        """
        outs = []
        for outp in outputs:
            Logger.log(f"Adding output '{outp.id()}' to '{self.identifier}'")
            node: OutputNode = OutputNode(outp)
            outs.append(node)
            self._add_node(node)
            self._outputs.append(node)
        return outs

    def _add_steps(self, steps: List[Step]) -> List[StepNode]:
        """
        Adds a list of steps to the graph and returns the created StepNodes
        """
        ss = []
        for step in steps:
            Logger.log(f"Adding step '{step.id()}' to '{self.identifier}'")
            node: StepNode = StepNode(step)

            self.has_subworkflow = self.has_subworkflow or step.tool().type() == ToolTypes.Workflow
            ss.append(node)
            self._add_node(node)
            self._steps.append(node)
        return ss

    def _add_node(self, node: Node):
        """
        Add a node to the graph, does the correct error checking
        """
        if node.id() in self._nodes:
            existing = self._nodes[node.id()]
            message = f"Attempted to add an input with the identifier '{node.id()}' but there was already " \
                f"an {existing.node_type} node with representation: {str(existing)}"
            Logger.log(message)
            raise DuplicateLabelIdentifier(message)
        self._nodes[node.id()] = node
        self.graph.add_node(node)

    @staticmethod
    def get_labels_and_node_by_inference(s) -> Tuple[str, Optional[Node]]:
        """
        Try to get the identifier and node of the
        :return: Tuple[str, Optional[Node]]
        """
        try:
            if type(s) == str:
                # just a string was passed,
                return s, None
            elif type(s) == tuple:
                # should have format (tag, Optional[Node]), check 105:step.py
                return s
            else:
                # s is hopefully a node (either InputNode | OutputNode)
                return s.id(), s
        except Exception as e:
            Logger.log_ex(e)
            return str(s), None

    def connect_inputs(self, step: Any, inputs: List[Any]) -> List[StepInput]:
        """
        Connect all the inputs to the step (infer the connections, will error if it can't be determined)
        :param step: Step to connect to
        :param inputs: Inputs to connect from
        :return: List[Edge]
        """
        return [self.add_piped_edge(step, i) for i in inputs]

    def add_edges(self, edges: List[Tuple[Any, Any]]):
        return [self.add_edge(e[0], e[1]) for e in edges]

    def add_pipe(self, *args) -> List[StepInput]:
        """
        Connect args[0] -> args[1] and args[1] to args[2] and ... and args[n-1] to args[n]
        :return:
        """
        if len(args) < 2:
            raise Exception("Must pass at least two properties to 'add_pipe'")

        # add_piped_edge removes the component on the start node (ie: n0 -> n1.tag, n1 -> n2.tag2, etc)
        return [self.add_piped_edge(args[i], args[i + 1]) for i in range(len(args) - 1)]

    def add_piped_edge(self, start, finish) -> StepInput:
        """
        add_piped_edge removes the component on the start node (ie: n0 -> n1.tag, n1 -> n2.tag2, etc)
        :param start: node to connect from
        :param finish: node to connect to
        :return: Edge between start and finish
        """
        s, sn = self.get_labels_and_node_by_inference(start)
        f, fn = self.get_labels_and_node_by_inference(finish)
        return self.add_edge((s.split("/")[0], sn), (f, fn))

    def add_default_value(self, tagged_node, default_value):
        """
        Adds a default value to the edge, will create an edge
        :param tagged_node: will be run through get_labels_and_node_by_inference: ( string | (string, node))
        :param default_value: value to be the default, should match the data_type
        """
        if not tagged_node:
            raise Exception("You must pass a non-optional node to 'add_default_value'")
        n, component = self.get_labels_and_node_by_inference(tagged_node)

        if not n:
            raise Exception(f"Couldn't resolve the node identifier from '{str(tagged_node)}'")
        components = n.split("/")
        node_id = components[0]
        node = self._nodes.get(node_id)
        if not node:
            raise Exception(f"Couldn't find a node in the graph with identifier '{node_id}', you must add the node "
                            f"to the graph to add a default value")
        if node.node_type != NodeTypes.TASK:
            raise Exception(f"You can only add a default value to a task (step) (not to '{node.id()}'")
        if not component:
            has_one_value = ", even if there's only one value" if len(node.inputs()) == 1 else ""
            raise Exception("You must fully qualify a node and it's components to add a default value" + has_one_value)

        tag = components[-1]

        # # Validate the type of data?
        # inp = node.inputs().get(tag)
        # inp.input_type == type(default)

        if tag in node.connection_map:
            e = node.connection_map[tag]
        else:
            e = StepInput(node, tag)
            node.connection_map[components[-1]] = e
        e.set_default(default_value)

    @staticmethod
    def validate_tag_parts(parts: List[str], tag: str, node_type: str):
        """
        Validates the parts in a tag (between 0 and 2)
        :param parts: list
        :param tag:
        :param node_type: start | finish for the error string
        :return:
        """
        if len(parts) == 0 or len(parts) > 2:
            message = f"The {node_type} tag '{tag}' contained an invalid number of parts ({len(parts)}), " \
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)

    def try_get_or_add_node(self, label: str, component: Optional[Node], node_type: str) -> Node:
        node = self._nodes.get(label)

        if node is not None:
            if component is not None and label != component.id():
                raise Exception(f"There already exists a node with id '{node.id()}' (clashes with {repr(component)}")
            return node

        Logger.log(f"Could not find {node_type} node with identifier '{label}' in the workflow")
        if component is not None:
            Logger.log(f"Adding '{component.id()}' to the workflow")
            return self._add_item(component)
        else:
            message = f"The node '{label}' was not found in the graph, please add the node to the graph " \
                f"and couldn't get it from the added component before trying to add an edge"
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)

    def add_edge(self, start, finish, default_value=None):
        """
        Add the start node as a source to the finish node, will create an edge if none exists
        :param start: an input or step node pair (via get_labels_and_node_by_inference)
        :param finish: a step or output node pair (via get_labels_and_node_by_inference)
        :param default_value: optionally provide a default value, if the start resolves to null
        :return: edge
        """

        s, s_component = self.get_labels_and_node_by_inference(start)
        f, f_component = self.get_labels_and_node_by_inference(finish)
        Logger.log(f"Adding edge between {s} and {f}")

        # set nodes
        s_parts = s.split("/")
        f_parts = f.split("/")
        stag, ftag = None, None

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
        self.validate_tag_parts(s_parts, s, "start")
        self.validate_tag_parts(f_parts, f, "finish")

        s_label = s_parts[0]
        f_label = f_parts[0]

        s_node = self.try_get_or_add_node(s_label, s_component, "start")
        f_node = self.try_get_or_add_node(f_label, f_component, "finish")

        s_node.set_depth(f_node.depth - 1)
        f_node.set_depth(s_node.depth + 1)

        if s_node.node_type == NodeTypes.INPUT:
            # Can guarantee the data_type
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(s_node, s_parts, referenced_by=f_parts)
            f_parts, f_type = self.get_tag_and_type_from_final_edge_node(f_node, f_parts, guess_type=s_type)
            stag, ftag = None, f_parts[-1]

        elif f_node.node_type == NodeTypes.OUTPUT:
            # We'll need to determine the data_type and assign it to the output
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(s_node, s_parts, referenced_by=f_parts)
            if s_type is None:
                keys = ", ".join(s_node.outputs().keys())
                raise Exception(f"The step '{s_node.id()}' must be fully qualified to connect to output node"
                                f" (expected one of: {keys})")

            step_has_scatter = any(e.has_scatter() for e in s_node.connection_map.values())

            if step_has_scatter:
                f_type = Array(s_type.received_type())
            else:
                f_type = s_type.received_type()

            f_node: OutputNode = f_node
            f_parts = [f_label]
            f_node.output.data_type = f_type
            stag, ftag = s_parts[-1], None
            Logger.log(f"Connecting '{s_node.id()}' to output '{f_node.id()}'")
        else:

            # We can't guarantee or enforce anything, so we'll see if we are fully qualified, and if not we can try
            # and guess by comparing types between the nodes, it doesn't matter which way we'll start though
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(s_node, s_parts, referenced_by=f_parts)
            f_parts, f_type = self.get_tag_and_type_from_final_edge_node(f_node, f_parts, guess_type=s_type)

            if s_type is None and f_type is not None:
                s_parts, s_type = self.get_tag_and_type_from_start_edge_node(
                    s_node,
                    s_parts,
                    referenced_by=f_parts,
                    guess_type=f_type
                )

            if s_type is None and f_type is None:
                # Try to guess both
                Logger.log(f"Trying to guess connection based on mutual types between "
                           f"'{s_node.id()}' and '{f_node.id()}'", LogLevel.WARNING)
                c = self.guess_connection_between_nodes(s_node, f_node)
                if c is not None:
                    (s_parts, s_type), (f_parts, f_type) = c

            stag, ftag = s_parts[-1], f_parts[-1]

        # We got one but not the other
        if s_type is not None and f_type is None:
            possible_outp_tags = ", ".join(f"{o.tag} ({o.input_type.name()})" for o in f_node.inputs().values())
            raise Exception(f"Couldn't definitively establish a connection from the start node '{s_node.id()}' "
                            f"with type {s_type.name()} to the finish node with possibilities: {possible_outp_tags}")
        if s_type is None and f_type is not None:
            possible_inp_tags = ", ".join(f"{i.tag} ({i.output_type.name()})" for i in s_node.outputs().values())
            raise Exception(f"Couldn't definitively establish a connection to the finish node '{f_node.id()}' "
                            f"with type {f_type.name()} from the start node with possibilities: {possible_inp_tags}")

        # We couldn't get either
        if s_type is None or f_type is None:
            ss = '/'.join(s_parts)
            ff = '/'.join(f_parts)
            raise Exception(f"Couldn't connect '{s_node.id()}' to '{f_node.id()}', failed to get types from "
                            f"tags (s: {ss} | f: {ff}) and couldn't uniquely isolate mutual types")

        # NOW: Let's build the connection (edge)

        if ftag in f_node.connection_map:
            step_inputs = f_node.connection_map[ftag]
        else:
            step_inputs = StepInput(f_node, ftag)
            f_node.connection_map[ftag] = step_inputs

        e = step_inputs.add_source(s_node, stag)
        if default_value:
            step_inputs.set_default(default_value)

        self.has_scatter = self.has_scatter or e.scatter
        self.has_multiple_inputs = self.has_multiple_inputs or step_inputs.multiple_inputs

        col = 'black' if e.compatible_types else 'r'
        self.graph.add_edge(s_node, f_node, type_match=e.compatible_types, color=col)
        return e

    def get_tag_and_type_from_start_edge_node(self, node: Node, input_parts: List[str], referenced_by: List[str],
                                              guess_type: Optional[DataType] = None) \
            -> Tuple[List[str], Optional[DataType]]:

        if guess_type is not None:
            guess_type = guess_type.received_type()

        if node.node_type == NodeTypes.OUTPUT:
            raise Exception(f"Can't join output node '{node.id()}' to '{referenced_by}'")

        lbl = input_parts[0]
        if node.node_type == NodeTypes.INPUT:
            if len(input_parts) != 1:
                s = "/".join(input_parts)
                message = f"The input tag '{s}' over-represents the INPUT node, this was corrected ({s} → {lbl})"
                Logger.log(message, LogLevel.WARNING)
                input_parts = [lbl]
            return input_parts, next(iter(node.outputs().values())).output_type

        # We are a step node now
        if not isinstance(node, StepNode):
            raise Exception("An internal error has occurred in 'get_tag_and_type_from_start_edge_node'"
                            " when converting Node to StepNode.")
        snode: StepNode = node
        outs: Dict[str, ToolOutput] = node.outputs()

        types = [(x, outs[x].output_type) for x in outs]

        if len(types) == 0:
            raise InvalidStepsException(f"The step '{referenced_by}' referenced the step '{snode.id()}' with tool "
                                        f"'{snode.step.tool().id()}' that has no outputs")
        elif len(types) == 1:
            tag = types[0][0]
            if len(input_parts) != 2:
                s = "/".join(input_parts)
                under_over = "under" if len(input_parts) < 2 else "over"
                Logger.warn(f"The node '{'/'.join(referenced_by)}' {under_over}-referenced an output of the tool "
                           f"'{snode.step.tool().id()}' (step: {node.id()}, this was automatically corrected "
                           f"({s} → {lbl}/{tag})")
            elif input_parts[-1] != tag:
                Logger.log(f"The node '{node.id()}' did not correctly reference an output of the tool "
                           f"'{snode.step.tool().id()}', this was automatically corrected "
                           f"({'/'.join(input_parts)} → {lbl}/{tag})", LogLevel.WARNING)

            input_parts = [lbl, tag]
            return input_parts, types[0][1]

        else:
            # if the edge tag doesn't match, we can give up
            if len(input_parts) < 2:
                if guess_type is not None:
                    # try to guess it from the outputs
                    compatible_types = [outs[x] for x in outs if guess_type.can_receive_from(outs[x].output_type)]
                    if len(compatible_types) == 1:
                        out = compatible_types[0]
                        Logger.info(
                            f"Guessed the compatible match for '{node.id()}' with source type '{out.output_type.id()}'"
                            f" → '{guess_type.id()}' ('{node.id()}/{out.tag}' → '{node.id()}')")
                        return [lbl, out.tag], out.output_type
                    else:
                        # Should this step also take into consideration the _optional_ nature of the node,
                        # ie: should it favour required nodes if that's an option
                        ultra_compatible = [outs[x] for x in outs if type(outs[x].output_type) == type(guess_type)]

                        if len(ultra_compatible) == 1:
                            ultra = ultra_compatible[0]
                            Logger.warn(f"There were {len(compatible_types)} matched types for the node '{node.id()}', "
                                        f"the program has guessed an exact compatible match of "
                                        f"type '{ultra.output_type.id()}' to tag '{ultra.tag}'")
                            return [lbl, ultra.tag], ultra.output_type
                        else:
                            s = "/".join(input_parts)
                            compat_str = ", ".join(f"{x.tag}: {x.output_type.id()}" for x in compatible_types)
                            raise Exception(
                                f"The node '{node.id()}' did not specify an input tag, and used '{guess_type.id()}'"
                                f" from the start node to guess the input by type, matching {len(compatible_types)}"
                                f" compatible ({compat_str}) and {len(ultra_compatible)} exact types."
                                f" You will need to provide more information to proceed.")
                else:
                    possible_tags = ", ".join(f"'{x}'" for x in outs)
                    s = "/".join(input_parts)
                    Logger.critical(f"The tag '{s}' could not uniquely identify an input of '{snode.id()}', requires the "
                                    f"one of the following tags: {possible_tags}")
                    return input_parts, None

            tag = input_parts[1]
            t = snode.step.tool().outputs_map().get(tag)
            if t:
                if len(input_parts) != 2:
                    Logger.log(f"The node '{node.id()}' did not correctly reference an output of the tool "
                               f"'{snode.step.tool().id()}', this was automatically corrected "
                               f"({'/'.join(input_parts)} → {lbl}/{tag})", LogLevel.WARNING)
                input_parts = [lbl, tag]
                return input_parts, t.output_type

            possible_tags = ", ".join(f"'{x}'" for x in outs)
            raise Exception(f"Could not identify an output called '{tag}' on the node '{node.id()}' "
                            f", possible tags: {possible_tags}")

        return input_parts, None

    @staticmethod
    def get_tag_and_type_from_final_edge_node(node: Node, input_parts: List[str],
                                              guess_type: Optional[DataType]) -> Tuple[List[str], Optional[DataType]]:
        if guess_type is not None:
            guess_type = guess_type.received_type()

        if node.node_type == NodeTypes.INPUT:
            raise Exception(f"Can't connect TO input node '{node.id()}'")
        if node.node_type == NodeTypes.OUTPUT:
            raise Exception("An internal error has occurred: output nodes should be filtered from the "
                            "function 'get_tag_and_type_from_final_edge_node'")

        if not isinstance(node, StepNode):
            raise Exception("An internal error has occurred in 'get_tag_and_type_from_start_edge_node'"
                            " when converting Node to StepNode.")

        snode: StepNode = node
        ins = snode.inputs()
        types = [(x, ins[x].input_type) for x in ins]
        lbl = input_parts[0]

        if len(types) == 0:
            raise InvalidStepsException(f"The step '{snode.id()}' has no inputs, and cannot be added to this workflow")

        if len(types) == 1:
            t = types[0]
            if len(input_parts) != 2:
                s = "/".join(input_parts)
                Logger.log(f"The node '{snode.id()}' did not correctly reference an input of the tool "
                           f"'{snode.step.tool().id()}', this was automatically corrected "
                           f"({s} → {lbl}/{t[0]})", LogLevel.WARNING)
            input_parts = [lbl, t[0]]
            return input_parts, t[1]
        elif len(input_parts) < 2:
            if guess_type is not None:
                # try to guess it from the outputs
                compatible_types = [ins[x] for x in ins if ins[x].input_type.can_receive_from(guess_type)]
                if len(compatible_types) == 1:
                    inp = compatible_types[0]
                    Logger.info(f"Guessed the compatible match for '{node.id()}' with source type '{guess_type.id()}'"
                                f" → '{inp.input_type.id()}' ('{node.id()}' → '{node.id()}/{inp.tag}')")
                    return [lbl, inp.tag], inp.input_type
                else:
                    # Should this step also take into consideration the _optional_ nature of the node,
                    # ie: should it favour required nodes if that's an option
                    ultra_compatible = [ins[x] for x in ins if type(ins[x].input_type) == type(guess_type)]

                    if len(ultra_compatible) == 1:
                        ultra = ultra_compatible[0]
                        Logger.warn(f"There were {len(compatible_types)} matched types for the node '{node.id()}', "
                                    f"the program has guessed an exact compatible match of "
                                    f"type '{ultra.input_type.id()}' to tag '{ultra.tag}'")
                        return [lbl, ultra.tag], ultra.input_type
                    else:
                        s = "/".join(input_parts)
                        compat_str = ", ".join(f"{x.tag}: {x.input_type.id()}" for x in compatible_types)
                        raise Exception(f"The end node '{node.id()}' did not specify where to join '{guess_type.id()}'"
                                        f" to, there were {len(compatible_types)} compatible types"
                                        f" ({compat_str}) and {len(ultra_compatible)} exact types."
                                        f" You will need to provide more information to proceed.")
            else:

                possible_tags = ", ".join(f"'{x}'" for x in ins)
                s = "/".join(input_parts)
                Logger.critical(f"The tag '{s}' could not uniquely identify an input of '{snode.id()}', requires the "
                                f"one of the following tags: {possible_tags}")
                return input_parts, None

        else:
            tag = input_parts[1]
            tool_input: ToolInput = snode.step.tool().inputs_map().get(tag)
            if tool_input:
                if len(input_parts) != 2:
                    s = "/".join(input_parts)
                    Logger.log(f"The node '{snode.id()}' did not correctly reference an input of the tool "
                               f"'{snode.step.tool().id()}', this was automatically corrected "
                               f"({s} → {lbl}/{tag})", LogLevel.WARNING)
                input_parts = [lbl, tag]
                return input_parts, tool_input.input_type

            possible_tags = ", ".join(f"'{x}'" for x in ins)
            raise Exception(f"Could not identify an input called '{tag}' on the node '{snode.id()}' "
                            f", possible tags: {possible_tags}")

        raise Exception("Unhandled pathway in 'get_tag_and_type_from_node'")

    @staticmethod
    def guess_connection_between_nodes(s_node: Node, f_node: Node) \
            -> Optional[Tuple[Tuple[List[str], DataType], Tuple[List[str], DataType]]]:
        outs, ins = s_node.outputs(), f_node.inputs()

        s_types: List[Tuple[List[str], DataType]] = [([s_node.id(), x], outs[x].output_type) for x in outs]
        f_types: List[Tuple[str, DataType]] = [([f_node.id(), x], ins[x].input_type) for x in ins]

        # O(n**2) for determining types
        matching_types: List[Tuple[Tuple[str, DataType], Tuple[str, DataType]]] = []
        for s_type in s_types:
            matching_types.extend([(s_type, f_type) for f_type in f_types if f_type[1].can_receive_from(s_type[1])])

        if len(matching_types) == 0:
            raise InvalidInputsException(f"Can't find a common type between '{s_node.id()}' and '{f_node.id()}'")

        if len(matching_types) > 1:
            compat = [str(t) for t in matching_types]
            message = f"Couldn't guess the single connection between '{s_node.id()}' and '{f_node.id()}', " \
                f"there was {len(compat)} compatible connections ({compat})"
            Logger.log(message, LogLevel.CRITICAL)
            return None

        matched = matching_types[0]
        Logger.info(f"Guessed the connection between nodes '{s_node.id()}")
        return matching_types[0]

    def cwl(self, is_nested_tool=False, with_docker=True):

        w = cwl.Workflow(self.identifier, self.label, self.doc)

        w.inputs = [i.cwl() for i in self._inputs]
        w.steps = [s.cwl(is_nested_tool=is_nested_tool) for s in self._steps]
        w.outputs = [o.cwl() for o in self._outputs]

        w.requirements.append(cwl.InlineJavascriptReq())

        if self.has_scatter:
            w.requirements.append(cwl.ScatterFeatureRequirement())
        if self.has_subworkflow:
            w.requirements.append(cwl.SubworkflowFeatureRequirement())
        if self.has_multiple_inputs:
            w.requirements.append(cwl.MultipleInputFeatureRequirement())

        tools = []
        tools_to_build: Dict[str, Tool] = {s.step.tool().id(): s.step.tool() for s in self._steps}
        for t in tools_to_build:
            tool: Tool = tools_to_build[t]
            if isinstance(tool, Workflow):
                wf_cwl, _, subtools = tool.cwl(is_nested_tool=True, with_docker=with_docker)
                tools.append(wf_cwl)
                tools.extend(subtools)
            else:
                tools.append(tool.cwl(with_docker=with_docker))

        inp = {i.id(): i.input.cwl_input() for i in self._inputs}

        return w.get_dict(), inp, tools

    def wdl2(self, with_docker=True):
        import wdlgen as wdl

        # Prepare for call aliases
        get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])

        tools: List[Tool] = [s.step.tool() for s in self._steps]
        tool_name_to_tool: Dict[str, Tool] = {t.id().lower(): t for t in tools}
        tool_name_to_alias = {}
        steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in self._steps}

        aliases = set()

        for tool in tool_name_to_tool:
            a = get_alias(tool).upper()
            s = a
            idx = 2
            while s in aliases:
                s = a + str(idx)
                idx += 1
            aliases.add(s)
            tool_name_to_alias[tool] = s

        # End call alias preparation

        w = wdl.Workflow(self.identifier)

        for i in self._inputs:
            wd = i.input.data_type.wdl2()
            if isinstance(wd, list):
                w.inputs.extend(wdl.Input(dt, i.id()) for dt in wd)
            else:
                w.inputs.append(wd)

        w.outputs = [wdl.Output(o.output.data_type.wdl2(), o.id(), "{a}.{b}".format(
            a=first_value(first_value(o.connection_map).source_map).start.id(),
            b=first_value(first_value(o.connection_map).source_map).stag
        )) for o in self._outputs]

        w.calls = [s.wdl2(s.step.tool().wdl_name(), s.id()) for s in self._steps]

        print(w.wdl())

    def wdl(self, **kwargs):

        get_alias = lambda t: t[0] + "".join([c for c in t[1:] if c.isupper()])

        tools: List[Tool] = [s.step.tool() for s in self._steps]
        tool_name_to_tool: Dict[str, Tool] = {t.id().lower(): t for t in tools}
        tool_name_to_alias = {}
        steps_to_alias: Dict[str, str] = {s.id().lower(): get_alias(s.id()).lower() for s in self._steps}

        aliases = set()

        for tool in tool_name_to_tool:
            a = get_alias(tool).upper()
            s = a
            idx = 2
            while s in aliases:
                s = a + str(idx)
                idx += 1
            aliases.add(s)
            tool_name_to_alias[tool] = s

        tab_char = '  '
        nline_char = '\n'

        import_str = "import \"tools/{tool_file}.wdl\" as {alias}"
        input_str = "{tb}{data_type} {identifier}"
        step_str = "{tb}call {tool_file}.{tool} as {alias} {{ input: {tool_mapping} }}"
        output_str = "{tb2}{data_type} {identifier} = {alias}.{outp}"

        imports = [import_str.format(
            tb=tab_char,
            tool_file=t,
            alias=tool_name_to_alias[t.lower()].upper()
        ) for t in tool_name_to_tool]

        inputs = [input_str.format(
            tb=tab_char,
            data_type=i.input.data_type.wdl(),
            identifier=i.id()
        ) for i in self._inputs]

        steps = [step_str.format(
            tb=tab_char,
            tool_file=tool_name_to_alias[s.step.tool().id().lower()].upper(),
            tool=s.step.tool().wdl_name(),
            alias=s.id(),
            tool_mapping=', '.join(s.wdl_map())  # [2 * tab_char + w for w in s.wdl_map()])
        ) for s in self._steps]

        outputs = [output_str.format(
            tb2=2 * tab_char,
            data_type=o.output.data_type.wdl(),
            identifier=o.id(),
            alias=first_value(first_value(o.connection_map).source_map).start.id(),
            outp=first_value(first_value(o.connection_map).source_map).stag
        ) for o in self._outputs]

        # imports = '\n'.join([f"import \"tools/{t}.wdl\" as {tool_name_to_alias[t.lower()].upper()}" for t in tool_name_to_tool])
        # inputs = '\n'.join([f"{tab_char}{i.input.data_type.wdl()} {i.id()}" for i in self._inputs])
        # steps = '\n'.join([f"{tab_char}call {tool_name_to_alias[s.step.get_tool().tool().lower()].upper()}"
        #                    f".{s.step.get_tool().wdl_name()} as {s.id()} {{input: \n{(',' + nline_char).join([2 * tab_char + w for w in s.wdl_map()])}\n{tab_char}}}" for s in self._steps])
        # outputs = '\n'.join([f"{2*tab_char}{o.output.data_type.wdl()} {o.id()} = {steps_to_alias[next(iter(o.connection_map.values()))[0].split('/')[0].lower()].lower()}.{o.id()}" for o in self._outputs])

        workflow = f"""
{nline_char.join(imports)}

workflow {self.identifier} {{
{nline_char.join(inputs)}

{nline_char.join(steps)}

{tab_char}output {{
{nline_char.join(outputs)}
{tab_char}}}
}}"""
        tools = {t.id(): t.wdl() for t in tools}

        inp = {f"{self.identifier}.{i.id()}": i.input.wdl_input() for i in self._inputs}

        return workflow, inp, tools

    def dump_cwl(self, to_disk: False, with_docker=True):
        import os, yaml, zipfile
        cwl_data, inp_data, tools_ar = self.cwl(with_docker=with_docker)

        d = os.path.expanduser("~") + f"/Desktop/{self.identifier}/cwl/"
        d_tools = d + "tools/"

        if not os.path.isdir(d):
            os.makedirs(d)
        if not os.path.isdir(d_tools):
            os.makedirs(d_tools)

        from cwlgen.cwlgen.utils import literal, literal_presenter
        yaml.add_representer(literal, literal_presenter)

        print(yaml.dump(cwl_data, default_flow_style=False))
        print(yaml.dump(inp_data, default_flow_style=False))
        for t in tools_ar:
            print(yaml.dump(t, default_flow_style=False))

        if to_disk:
            os.chdir(d)
            with open(d + self.id() + ".cwl", "w+") as cwl:
                Logger.log(f"Writing {self.identifier}.cwl to disk")
                yaml.dump(cwl_data, cwl, default_flow_style=False)
                Logger.log(f"Written {self.identifier}.cwl to disk")

            with open(d + self.id()+ "-job.yml", "w+") as cwl:
                Logger.log(f"Writing {self.identifier}-job.yml to disk")
                yaml.dump(inp_data, cwl, default_flow_style=False)
                Logger.log(f"Written {self.identifier}-job.yml to disk")

            # z = zipfile.ZipFile(d + "tools.zip", "w")
            for tool in tools_ar:
                tool_name = tool["id"]
                tool_filename = tool_name + ".cwl"
                with open(d_tools + tool_filename, "w+") as cwl:
                    Logger.log(f"Writing {tool_name}.cwl to disk")
                    yaml.dump(tool, cwl, default_flow_style=False)
                    Logger.log(f"Written {tool_name}.cwl to disk")

            #     z.write("tools/" + tool_filename)
            # z.close()



    def dump_wdl(self, to_disk: False):
        import os, json
        _ = self.wdl2()
        return
        wdl_data, inp_data, tools_dict = self.wdl2()
        print(wdl_data)
        print("================")
        print(inp_data)
        print("================")
        print("\n*******\n".join(tools_dict.values()))

        d = os.path.expanduser("~") + f"/Desktop/{self.identifier}/wdl/"
        d_tools = d + "tools/"
        if not os.path.isdir(d):
            os.makedirs(d)
        if not os.path.isdir(d_tools):
            os.makedirs(d_tools)

        if to_disk:
            with open(d + self.identifier + ".wdl", "w+") as wdl:
                Logger.log(f"Writing {self.identifier}.wdl to disk")
                wdl.write(wdl_data)
                Logger.log(f"Written {self.identifier}.wdl to disk")

            with open(d + self.identifier + "-job.json", "w+") as inp:
                Logger.log(f"Writing {self.identifier}-job.json to disk")
                json.dump(inp_data, inp)
                Logger.log(f"Written {self.identifier}-job.json to disk")

            for tool_name in tools_dict:
                tool = tools_dict[tool_name]
                with open(d_tools + tool_name + ".wdl", "w+") as wdl:
                    Logger.log(f"Writing {tool_name}.cwl to disk")
                    wdl.write(tool)
                    Logger.log(f"Written {tool_name}.cwl to disk")
