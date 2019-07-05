from typing import Dict, List, Tuple, Optional, Any, Union, Iterable

import os
import networkx as nx

import janis.translations as translations
from janis.graph.node import Node, NodeTypes, layout_nodes2
from janis.graph.stepinput import StepInput
from janis.tool.commandtool import CommandTool
from janis.tool.tool import Tool, ToolInput, ToolOutput, ToolTypes, ToolType
from janis.translations import ExportPathKeywords
from janis.types.common_data_types import Array
from janis.types.data_types import DataType
from janis.utils.errors import (
    DuplicateLabelIdentifier,
    InvalidNodeIdentifier,
    NodeNotFound,
    InvalidStepsException,
    InvalidInputsException,
)
from janis.utils.logger import Logger, LogLevel
from janis.utils.metadata import WorkflowMetadata
from janis.utils.validators import Validators
from janis.workflow.input import Input, InputNode
from janis.workflow.output import Output, OutputNode
from janis.workflow.step import Step, StepNode


class Workflow(Tool):
    """
    The Workflow contains information about inputs, outputs, steps and their relations to each other.

    You can construct workflows in 2 ways:

        1. Create an instance of Workflow

        2. Subclass Workflow (if you're packaging a workflow as a tool, this is the preferred approach).


    """

    def __init__(
        self, identifier: str, friendly_name: str = None, doc: Optional[str] = None
    ):
        """
        Initialise the workflow
        :param identifier: uniquely identifies the workflow
        :param friendly_name: a label that the engine may use to represent the workflow / used as name in docs
        :param doc: Documentation of the workflow, what is
        """
        Logger.log(f"Creating workflow with identifier: '{identifier}'")

        if not Validators.validate_identifier(identifier):
            raise Exception(
                f"The identifier '{identifier}' was not validated by '{Validators.identifier_regex}' "
                f"(must start with letters, and then only contain letters, numbers and an underscore)"
            )

        self.identifier = identifier
        self.name = friendly_name
        self._metadata = WorkflowMetadata(documentation=doc)

        # The following variables allow us to quickly check data about the graph

        self._nodes: Dict[str, Node] = {}  # Look up a node by its identifier

        self._inputs: List[InputNode] = []
        self._steps: List[StepNode] = []
        self._outputs: List[OutputNode] = []

        # We store approximate relationships on the MultiDiGraph, however
        # it can't store multiple connections between the same pairs of nodes very easily
        self.graph: nx.MultiDiGraph = nx.MultiDiGraph()

        # Flags for different requirements that a workflow might need
        self.has_scatter = False
        self.has_subworkflow = False
        self.has_multiple_inputs = False

    def id(self) -> str:
        """
        Returns the identifier of the workflow
        """
        return self.identifier

    def friendly_name(self):
        return self.name

    def metadata(self):
        return self._metadata

    def doc(self):
        return self._metadata.documentation

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
        return [
            ToolOutput(o.id(), o.output.data_type, o.output.doc)
            for o in self._outputs
            if o.output.data_type is not None
        ]

    # GRAPH CONSTRUCTION

    def draw_graph(self):
        """
        Use matplotlib to draw the a diagram of the workflow (currently neglects subworkflows).
        This will pause the execution of the program until the graph is closed.
        """
        import matplotlib.pyplot as plt

        default_color = "black"

        G = self.graph
        edges_attributes = [G.edges[e] for e in G.edges]
        edge_colors = [
            x["color"] if "color" in x else default_color for x in edges_attributes
        ]
        node_colors = [NodeTypes.to_col(x.node_type) for x in G.nodes]

        pos = layout_nodes2(
            list(G.nodes)
        )  # Manual layout engine of the graph, need to make this better

        for n in pos:
            # noinspection PyUnresolvedReferences
            G.node[n]["pos"] = pos[n]

        nx.draw(
            G, pos=pos, edge_color=edge_colors, node_color=node_colors, with_labels=True
        )
        plt.show()

    # For the most part, Just add edges as this will automatically add nodes to the graph

    def add_items(self, *args) -> List[Any]:
        """
        Supports adding (Input | Step | Output | Tuple = construct edge)
        """
        mixed: Iterable[Any] = args

        if isinstance(args[0], list):
            if len(args) > 1:
                raise Exception("Invalid type of arguments")
            mixed = args[0]

        return [self._add_item(i) for i in mixed]

    def _add_item(self, item: Union[Input, Step, Output, Tuple]):
        """
        Generic add whatever into the workflow: (Input | Step | Output | Tuple = construct edge)
        :return:
        """
        if isinstance(item, Input):
            return self._add_input(item)
        elif isinstance(item, Step):
            return self._add_steps(item)
        elif isinstance(item, Output):
            return self._add_outputs(item)
        elif isinstance(item, tuple):
            return self.add_edge(item[0], item[1])
        else:
            raise Exception(f"Unexpected type '{type(item)}' passed to 'add_nodes'")

    def _add_input(self, inp: Input) -> InputNode:
        """
        Add an input into the graph, and returns the created InputNode
        :param: inp
        :type: janis.Input
        """
        Logger.log(f"Adding input '{inp.id()}' to '{self.identifier}'")
        node: InputNode = InputNode(inp)
        self._add_node(node)
        self._inputs.append(node)
        return node

    def _add_outputs(self, outp: Output) -> OutputNode:
        """
        Adds an output to the graph, and returns the created OutputNode
        :param: outp
        :type: janis.Output
        """
        Logger.log(f"Adding output '{outp.id()}' to '{self.identifier}'")
        node: OutputNode = OutputNode(outp)
        self._add_node(node)
        self._outputs.append(node)
        return node

    def _add_steps(self, step: Step) -> StepNode:
        """
        Adds a list of steps to the graph and returns the created StepNodes
        :param: step
        :type: janis.Step
        """
        Logger.log(f"Adding step '{step.id()}' to '{self.identifier}'")
        node: StepNode = StepNode(step)

        self.has_subworkflow = (
            self.has_subworkflow or step.tool().type() == ToolTypes.Workflow
        )
        self._add_node(node)
        self._steps.append(node)
        return node

    def _add_node(self, node: Node) -> Node:
        """
        Add a node to the graph. This won't allow you to add the same node twice
        :param: node
        :type: janis.Node
        """
        if node.id() in self._nodes:
            existing = self._nodes[node.id()]
            message = (
                f"Attempted to add an input with the identifier '{node.id()}' but there was already "
                f"an {existing.node_type} node with representation: {str(existing)}"
            )
            Logger.log(message)
            raise DuplicateLabelIdentifier(message)
        self._nodes[node.id()] = node
        self.graph.add_node(node)
        return node

    @staticmethod
    def get_labels_and_node_by_inference(s) -> Tuple[str, Optional[Node]]:
        """
        Try to get the identifier and node of what's passed to the add_edges method
        :return: Tuple[str, Optional[Node]]
        """
        try:
            if type(s) == str:
                # just a string was passed,
                return s, None
            elif type(s) == tuple:
                # should have format (tag, Optional[Node]), check 105:step.py
                return s[0], s[1]
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
            message = (
                f"The {node_type} tag '{tag}' contained an invalid number of parts ({len(parts)}), "
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            )
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)

    def try_get_or_add_node(
        self,
        identifier: str,
        component: Optional[Union[Input, Step, Output]],
        node_type: str,
    ) -> Node:
        """
        Given the label (and the component [Input \ Step | Output], get the node if it's in the graph
        (and check that it's referencing the correct component), else add it to the graph.
        :param identifier: identifier
        :type identifier: str
        :param component: the step | input | output that the node contains
        :type component: Input | Step | Output
        :param node_type: Used to better classify errors ["input" | "step" | "output"]
        :type node_type: str
        :return: Node | Exception
        :returns: Node
        """
        node = self._nodes.get(identifier)

        if node is not None:
            # Try to do some proper matching
            node_comp: Optional[Any] = None
            if isinstance(node, InputNode):
                node_comp = node.input
            elif isinstance(node, StepNode):
                node_comp = node.step
            elif isinstance(node, OutputNode):
                node_comp = node.output

            if component is not None and (
                identifier != component.id()
                or (
                    isinstance(component, Input)
                    and not node.node_type == NodeTypes.INPUT
                )
                or (
                    isinstance(component, Step) and not node.node_type == NodeTypes.TASK
                )
                or (
                    isinstance(component, Output)
                    and not node.node_type == NodeTypes.OUTPUT
                )
                or (node_comp != component)
            ):
                raise Exception(
                    f"There already exists a node (and component) with id '{node.id()}'. The added "
                    f"component ('{repr(component)}') clashes with '{repr(node_comp)}')."
                )
            return node

        Logger.log(
            f"Could't find a(n) {node_type} node with the identifier '{identifier}' in the workflow."
        )
        if component is not None:
            Logger.log(f"Adding '{component.id()}' to the workflow")
            return self._add_item(component)
        else:
            message = (
                f"There was no node or component referenced by '{identifier}' in the graph. When creating edges "
                "by identifiers, you must add the component to the graph first"
            )
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)

    def add_edge(self, start, finish):
        """
        Add the start node as a source to the finish node, will create an edge if none exists
        :param start: an input or step node pair (via get_labels_and_node_by_inference)
        :param finish: a step or output node pair (via get_labels_and_node_by_inference)
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
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(
                s_node, s_parts, referenced_by=f_parts
            )
            f_parts, f_type = self.get_tag_and_type_from_final_edge_node(
                f_node, f_parts, referenced_from=s_parts, guess_type=s_type
            )
            stag, ftag = None, f_parts[-1]

        elif f_node.node_type == NodeTypes.OUTPUT:
            # We'll need to determine the data_type and assign it to the output
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(
                s_node, s_parts, referenced_by=f_parts
            )
            if s_type is None:
                keys = ", ".join(s_node.outputs().keys())
                raise Exception(
                    f"The step '{s_node.id()}' must be fully qualified to connect to output node"
                    f" (expected one of: {keys})"
                )

            step_has_scatter = any(
                e.has_scatter() for e in s_node.connection_map.values()
            )

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
            s_parts, s_type = self.get_tag_and_type_from_start_edge_node(
                s_node, s_parts, referenced_by=f_parts
            )
            f_parts, f_type = self.get_tag_and_type_from_final_edge_node(
                f_node, f_parts, referenced_from=s_parts, guess_type=s_type
            )

            if s_type is None and f_type is not None:
                s_parts, s_type = self.get_tag_and_type_from_start_edge_node(
                    s_node, s_parts, referenced_by=f_parts, guess_type=f_type
                )

            if s_type is None and f_type is None:
                # Try to guess both
                Logger.info(
                    f"Trying to guess connection based on mutual types between "
                    f"'{s_node.id()}' and '{f_node.id()}'"
                )
                c = self.guess_connection_between_nodes(s_node, f_node)
                if c is not None:
                    (s_parts, s_type), (f_parts, f_type) = c

            stag, ftag = s_parts[-1], f_parts[-1]

        # We got one but not the other
        if s_type is not None and f_type is None:
            possible_outp_tags = ", ".join(
                f"{o.tag} ({o.input_type.name()})" for o in f_node.inputs().values()
            )
            raise Exception(
                f"Couldn't definitively establish a connection from the start node '{s_node.id()}' "
                f"with type {s_type.name()} to the finish node with possibilities: {possible_outp_tags}"
            )
        if s_type is None and f_type is not None:
            possible_inp_tags = ", ".join(
                f"{i.tag} ({i.output_type.name()})" for i in s_node.outputs().values()
            )
            raise Exception(
                f"Couldn't definitively establish a connection to the finish node '{f_node.id()}' "
                f"with type {f_type.name()} from the start node with possibilities: {possible_inp_tags}"
            )

        # We couldn't get either
        if s_type is None or f_type is None:
            ss = "/".join(s_parts)
            ff = "/".join(f_parts)
            raise Exception(
                f"Couldn't connect '{s_node.id()}' to '{f_node.id()}', failed to get types from "
                f"tags (s: {ss} | f: {ff}) and couldn't uniquely isolate mutual types"
            )

        # NOW: Let's build the connection (edge)

        if ftag in f_node.connection_map:
            step_inputs = f_node.connection_map[ftag]
        else:
            step_inputs = StepInput(f_node, ftag)
            f_node.connection_map[ftag] = step_inputs

        e = step_inputs.add_source(s_node, stag)

        self.has_scatter = self.has_scatter or e.scatter
        self.has_multiple_inputs = (
            self.has_multiple_inputs or step_inputs.multiple_inputs
        )

        col = "black" if e.compatible_types else "r"
        self.graph.add_edge(s_node, f_node, type_match=e.compatible_types, color=col)
        return e

    def get_tag_and_type_from_start_edge_node(
        self,
        node: Node,
        input_parts: List[str],
        referenced_by: List[str],
        guess_type: Optional[DataType] = None,
    ) -> Tuple[List[str], Optional[DataType]]:

        if guess_type is not None:
            guess_type = guess_type.received_type()

        if node.node_type == NodeTypes.OUTPUT:
            ref_by = ".".join(referenced_by)
            raise Exception(
                f"Joining an INPUT node ('{node.id()}') to an  OUTPUT node ('{ref_by}') "
                f"has been explicitly disabled"
            )

        lbl = input_parts[0]
        if node.node_type == NodeTypes.INPUT:
            if len(input_parts) != 1:
                # You can't really get here if you're building pipelines with the Python classes, but you can build
                # the graph purely by referencing strings, so this is still needed.

                s = ".".join(input_parts)
                f = ".".join(referenced_by)
                message = f"The step '{f}' over-referenced the INPUT node, this was corrected ({s} → {lbl})"
                Logger.log(message, LogLevel.WARNING)
                input_parts = [lbl]
            return input_parts, next(iter(node.outputs().values())).output_type

        # We are a step node now
        if not isinstance(node, StepNode):
            raise Exception(
                "An internal error has occurred in 'get_tag_and_type_from_start_edge_node'"
                " when converting Node to StepNode."
            )
        snode: StepNode = node
        outs: Dict[str, ToolOutput] = node.outputs()

        types = [(x, outs[x].output_type) for x in outs]

        if len(types) == 0:
            raise InvalidStepsException(
                f"The step '{referenced_by}' referenced the step '{snode.id()}' with tool "
                f"'{snode.step.tool().id()}' that has no outputs"
            )
        elif len(types) == 1:
            tag = types[0][0]
            if len(input_parts) != 2:
                s = "/".join(input_parts)
                under_over = "under" if len(input_parts) < 2 else "over"
                Logger.info(
                    f"The node '{'.'.join(referenced_by)}' {under_over}-referenced an output the step "
                    f"'{node.id()}' (tool: '{snode.step.tool().id()}'), this was automatically corrected "
                    f"({s} → {lbl}.{tag})"
                )

            elif input_parts[-1] != tag:
                Logger.critical(
                    f"The node '{node.id()}' did not correctly reference an output of the step "
                    f"'{node.id()}' (tool: '{snode.step.tool().id()}'), this was automatically corrected "
                    f"({'/'.join(input_parts)} → {lbl}/{tag})"
                )

            input_parts = [lbl, tag]
            return input_parts, types[0][1]

        else:
            # if the edge tag doesn't match, we can give up
            if len(input_parts) < 2:
                if guess_type is not None:
                    # try to guess it from the outputs
                    compatible_types = [
                        outs[x]
                        for x in outs
                        if (not outs[x].output_type.optional)
                        and guess_type.can_receive_from(outs[x].output_type)
                    ]
                    if len(compatible_types) == 1:
                        out = compatible_types[0]
                        Logger.info(
                            f"Guessed the compatible match for '{node.id()}' with source type '{out.output_type.id()}'"
                            f" → '{guess_type.id()}' ('{node.id()}/{out.tag}' → '{node.id()}')"
                        )
                        return [lbl, out.tag], out.output_type
                    else:
                        # Should this step also take into consideration the _optional_ nature of the node,
                        # ie: should it favour required nodes if that's an option
                        ultra_compatible = [
                            outs[x]
                            for x in outs
                            if (not outs[x].output_type.optional)
                            and type(outs[x].output_type) == type(guess_type)
                        ]

                        if len(ultra_compatible) == 1:
                            ultra = ultra_compatible[0]
                            Logger.warn(
                                f"There were {len(compatible_types)} matched types for the node '{node.id()}', "
                                f"the program has guessed an exact compatible match of "
                                f"type '{ultra.output_type.id()}' to tag '{ultra.tag}'"
                            )
                            return [lbl, ultra.tag], ultra.output_type
                        else:
                            s = "/".join(input_parts)
                            compat_str = ", ".join(
                                f"{x.tag}: {x.output_type.id()}"
                                for x in compatible_types
                            )
                            raise Exception(
                                f"The node '{node.id()}' did not specify an input tag, and used '{guess_type.id()}'"
                                f" from the start node to guess the input by type, matching {len(compatible_types)}"
                                f" compatible ({compat_str}) and {len(ultra_compatible)} exact types."
                                f" You will need to provide more information to proceed."
                            )
                else:
                    possible_tags = ", ".join(f"'{x}'" for x in outs)
                    s = "/".join(input_parts)
                    Logger.critical(
                        f"The tag '{s}' could not uniquely identify an input of '{snode.id()}', requires the "
                        f"one of the following tags: {possible_tags}"
                    )
                    return input_parts, None

            tag = input_parts[1]
            t = snode.step.tool().outputs_map().get(tag)
            if t:
                if len(input_parts) != 2:
                    Logger.log(
                        f"The node '{node.id()}' did not correctly reference an output of the tool "
                        f"'{snode.step.tool().id()}', this was automatically corrected "
                        f"({'/'.join(input_parts)} → {lbl}/{tag})",
                        LogLevel.WARNING,
                    )
                input_parts = [lbl, tag]
                return input_parts, t.output_type

            possible_tags = ", ".join(f"'{x}'" for x in outs)
            raise Exception(
                f"Could not identify an output called '{tag}' on the node '{node.id()}' "
                f", possible tags: {possible_tags}"
            )

    @staticmethod
    def get_tag_and_type_from_final_edge_node(
        node: Node,
        input_parts: List[str],
        guess_type: Optional[DataType],
        referenced_from: List[str],
    ) -> Tuple[List[str], Optional[DataType]]:
        if guess_type is not None:
            guess_type = guess_type.received_type()

        if node.node_type == NodeTypes.INPUT:
            raise Exception(f"Can't connect TO input node '{node.id()}'")
        if node.node_type == NodeTypes.OUTPUT:
            raise Exception(
                "An internal error has occurred: output nodes should be filtered from the "
                "function 'get_tag_and_type_from_final_edge_node'"
            )

        if not isinstance(node, StepNode):
            raise Exception(
                "An internal error has occurred in 'get_tag_and_type_from_start_edge_node'"
                " when converting Node to StepNode."
            )

        snode: StepNode = node
        ins = snode.inputs()
        types = [(x, ins[x].input_type) for x in ins]
        lbl = input_parts[0]

        if len(types) == 0:
            raise InvalidStepsException(
                f"The step '{snode.id()}' has no inputs, and cannot be added to this workflow"
            )

        if len(types) == 1:
            t = types[0]
            if len(input_parts) != 2:
                s = "/".join(input_parts)
                Logger.info(
                    f"The node '{snode.id()}' was not a fully qualified input of the tool "
                    f"'{snode.step.tool().id()}', this was automatically corrected "
                    f"({s} → {lbl}.{t[0]})"
                )
            input_parts = [lbl, t[0]]
            return input_parts, t[1]
        elif len(input_parts) < 2:
            if guess_type is not None:
                # try to guess it from the outputs
                compatible_types = [
                    ins[x]
                    for x in ins
                    if (not ins[x].input_type.optional)
                    and ins[x].input_type.can_receive_from(guess_type)
                ]
                if len(compatible_types) == 1:
                    inp = compatible_types[0]
                    Logger.info(
                        f"Guessed the compatible match between {'.'.join(referenced_from)} → {node.id()}.[unknown] with source type '{guess_type.id()}'"
                        f" → '{inp.input_type.id()}' ('{node.id()}' → '{node.id()}.{inp.tag}')"
                    )
                    return [lbl, inp.tag], inp.input_type
                else:
                    # Should this step also take into consideration the _optional_ nature of the node,
                    # ie: should it favour required nodes if that's an option
                    ultra_compatible = [
                        ins[x]
                        for x in ins
                        if (not ins[x].input_type.optional)
                        and type(ins[x].input_type) == type(guess_type)
                    ]

                    if len(ultra_compatible) == 1:
                        ultra = ultra_compatible[0]
                        Logger.warn(
                            f"There were {len(compatible_types)} matched types for the node '{node.id()}', "
                            f"the program has guessed an exact compatible match of "
                            f"type '{ultra.input_type.id()}' to tag '{ultra.tag}'"
                        )
                        return [lbl, ultra.tag], ultra.input_type
                    else:
                        s = ".".join(input_parts)
                        compat_str = ", ".join(
                            f"{x.tag}: {x.input_type.id()}" for x in compatible_types
                        )
                        raise Exception(
                            f"An error occurred when connecting '{'.'.join(referenced_from)}' to "
                            f"'{'.'.join(input_parts)}'. There were {len(compatible_types)} compatible "
                            f"types ({compat_str}) and {len(ultra_compatible)} exact types."
                            f" You will need to provide more information to proceed."
                        )
            else:

                possible_tags = ", ".join(f"'{x}'" for x in ins)
                s = ".".join(input_parts)
                Logger.critical(
                    f"The tag '{s}' could not uniquely identify an input of '{snode.id()}', requires the "
                    f"one of the following tags: {possible_tags}"
                )
                return input_parts, None

        else:
            tag = input_parts[1]
            tool_input: Optional[ToolInput] = snode.step.tool().inputs_map().get(tag)
            if tool_input:
                if len(input_parts) != 2:
                    s = "/".join(input_parts)
                    Logger.log(
                        f"The node '{snode.id()}' did not correctly reference an input of the tool "
                        f"'{snode.step.tool().id()}', this was automatically corrected "
                        f"({s} → {lbl}/{tag})",
                        LogLevel.WARNING,
                    )
                input_parts = [lbl, tag]
                return input_parts, tool_input.input_type

            possible_tags = ", ".join(f"'{x}'" for x in ins)
            raise Exception(
                f"Could not identify an input called '{tag}' on the node '{snode.id()}' "
                f", possible tags: {possible_tags}"
            )

        # raise Exception("Unhandled pathway in 'get_tag_and_type_from_node'")

    @staticmethod
    def guess_connection_between_nodes(
        s_node: Node, f_node: Node
    ) -> Optional[Tuple[Tuple[List[str], DataType], Tuple[List[str], DataType]]]:
        outs, ins = s_node.outputs(), f_node.inputs()

        s_types: List[Tuple[List[str], DataType]] = [
            ([s_node.id(), x], outs[x].output_type) for x in outs
        ]
        f_types: List[Tuple[List[str], DataType]] = [
            ([f_node.id(), x], ins[x].input_type) for x in ins
        ]

        # O(n**2) for determining types
        matching_types: List[
            Tuple[Tuple[List[str], DataType], Tuple[List[str], DataType]]
        ] = []
        for s_type in s_types:
            matching_types.extend(
                [
                    (s_type, f_type)
                    for f_type in f_types
                    if f_type[1].can_receive_from(s_type[1])
                ]
            )

        if len(matching_types) == 0:
            raise InvalidInputsException(
                f"Can't find a common type between '{s_node.id()}' and '{f_node.id()}'"
            )

        if len(matching_types) > 1:
            compat = [str(t) for t in matching_types]
            message = (
                f"Couldn't guess the single connection between '{s_node.id()}' and '{f_node.id()}', "
                f"there was {len(compat)} compatible connections ({compat})"
            )
            Logger.log(message, LogLevel.CRITICAL)
            return None

        matched = matching_types[0]
        Logger.info(f"Guessed the connection between nodes '{s_node.id()}")
        return matched

    # TRANSLATIONS

    def translate(
        self,
        translation: translations.SupportedTranslation,
        to_console=True,
        to_disk=False,
        write_inputs_file=True,
        with_docker=True,
        with_hints=False,
        with_resource_overrides=False,
        should_validate=False,
        should_zip=True,
        export_path=ExportPathKeywords.default,
        merge_resources=False,
        hints=None,
        allow_null_if_not_optional=True,
        additional_inputs: Dict = None,
    ):
        return translations.translate_workflow(
            self,
            translation=translation,
            to_console=to_console,
            to_disk=to_disk,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
            should_zip=should_zip,
            export_path=export_path,
            write_inputs_file=write_inputs_file,
            should_validate=should_validate,
            merge_resources=merge_resources,
            hints=hints,
            allow_null_if_not_optional=allow_null_if_not_optional,
            additional_inputs=additional_inputs,
        )

    def generate_resources_file(
        self,
        translation: translations.SupportedTranslation,
        hints: Dict[str, Any] = None,
        to_console=True,
    ):
        tr = translations.build_resources_input(self, translation, hints)
        if to_console:
            print(tr)
        return tr

    def get_tools(self) -> Dict[str, CommandTool]:
        tools: Dict[str, CommandTool] = {}
        for t in self._steps:
            tl = t.step.tool()
            if isinstance(tl, Workflow):
                tools.update(tl.get_tools())
            elif t.id() not in tools:
                tools[tl.id()] = tl
        return tools

    def generate_resources_table(
        self,
        hints: Dict[str, Any],
        to_console=True,
        to_disk=False,
        output_type: str = "tsv",
    ):
        delim = "\t" if output_type == "tsv" else ","

        tools = self.get_tools()
        header = ["name", "cpu", "memory (GB)"]
        data = []

        for t in sorted(tools.keys()):
            tool = tools[t]
            data.append([tool.id(), tool.cpus(hints), tool.memory(hints)])

        data.sort(key=lambda a: a[0].lower())

        data.insert(0, header)

        if to_console:
            import tabulate

            print(tabulate.tabulate(data, headers="firstrow"))
        if to_disk:
            import csv

            d = ExportPathKeywords.resolve(
                ExportPathKeywords.default_no_spec, None, self.id()
            )
            path = d + f"resources.{output_type}"

            if not os.path.isdir(d):
                os.makedirs(d)

            Logger.info(f"Writing resources {output_type} to '{path}'")
            with open(path, "w+") as mf:
                writer = csv.writer(mf, delimiter=delim)
                for row in data:
                    writer.writerow(row)

        return data
