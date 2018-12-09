from typing import Dict, List, Tuple, Set

from Pipeline.types.data_types import DataType
from Pipeline.utils.logger import Logger, LogLevel
from Pipeline.utils.descriptor import BaseDescriptor, GetOnlyDescriptor

import networkx as nx

from Pipeline.graph.node import Node, NodeTypes, NodeAndTag
from Pipeline.utils.errors import DuplicateLabelIdentifier, InvalidNodeIdentifier, NodeNotFound, InvalidStepsException, \
    InvalidInputsException
from Pipeline.workflow.input import Input, InputNode
from Pipeline.workflow.output import Output, OutputNode
from Pipeline.workflow.step import Step, StepNode
from Pipeline.tool.tool import Tool, ToolInput, ToolOutput


class Workflow:

    name = BaseDescriptor()

    def __init__(self, name: str):
        Logger.log(f"Creating workflow with name: '{name}'")

        self.name = name

        self._labels: Dict[str, Node] = {}

        self._inputs: List[InputNode] = []
        self._steps: List[StepNode] = []
        self._outputs: List[OutputNode] = []

        self.graph: nx.MultiDiGraph = nx.MultiDiGraph()
        self.connections: Dict[str, Set[NodeAndTag]] = {}

    def add_input(self, inp: Input):
        Logger.log(f"Adding input '{inp.id()}' to '{self.name}'")
        node: InputNode = InputNode(inp)
        self._add_node(node)
        self._inputs.append(node)

    def add_output(self, outp: Output):
        Logger.log(f"Adding output '{outp.id()}' to '{self.name}'")
        node: OutputNode = OutputNode(outp)
        self._add_node(node)
        self._outputs.append(node)

    def add_step(self, step: Step):
        Logger.log(f"Adding step '{step.id()}' to '{self.name}'")
        node: StepNode = StepNode(step)
        self._add_node(node)
        self._steps.append(node)

    def _add_node(self, node: Node):

        if node.id() in self._labels:
            existing = self._labels[node.id()]
            message = f"Attempted to add an input with the identifier '{node.id()}' but there was already " \
                f"an {existing.node_type} node with representation: {str(existing)}"
            Logger.log(message)
            raise DuplicateLabelIdentifier(message)
        self._labels[node.id()] = node
        self.graph.add_node(node)

    def _remove_node(self, node: Node):

        if node.id() not in self._labels:
            Logger.log(f"Couldn't find node with id '{node.id()}' to remove")
            return

    @staticmethod
    def get_label_by_inference(s) -> str:
        try:
            if type(s) == str:
                return s
            else:
                return s.id()
        except Exception as e:
            Logger.log(e, LogLevel.CRITICAL)

    def add_edge(self, start, finish):

        s = self.get_label_by_inference(start)
        f = self.get_label_by_inference(finish)
        Logger.log(f"Adding edge between {s} and {f}", LogLevel.INFO)

        # set nodes
        s_parts = s.split("/")
        f_parts = f.split("/")

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

        if len(s_parts) == 0 or len(s_parts) > 2:
            message = f"The start tag '{s}' contained an invalid number of parts ({len(s_parts)}), " \
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)
        if len(f_parts) == 0 or len(f_parts) > 2:
            message = f"The finish tag '{f}' contained an invalid number of parts ({len(f_parts)}), " \
                f"it should follow the format: 'label/tag' (1 or 2 parts)"
            Logger.log(message, LogLevel.CRITICAL)
            raise InvalidNodeIdentifier(message)

        s_label = s_parts[0]
        f_label = f_parts[0]

        s_node = self._labels.get(s_label)
        f_node = self._labels.get(f_label)

        if s_node is None:
            message = f"The node '{s_label}' was not found in the graph, " \
                f"please add the node to the graph before trying to add an edge"
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)
        if f_node is None:
            message = f"The node '{f_label}' was not found in the graph, " \
                f"please add the node to the graph before trying to add an edge"
            Logger.log(message, LogLevel.CRITICAL)
            raise NodeNotFound(message)

        # we have the two nodes

        s_node.set_depth(f_node.depth - 1)
        f_node.set_depth(s_node.depth + 1)

        s_type = self.get_tag_and_type_from_node(s_node, s_parts, True)

        if f_node.node_type == NodeTypes.OUTPUT:
            f_type = f_parts[0], s_type[1]
            Logger.log("Connecting to output")
        else:
            f_type = self.get_tag_and_type_from_node(f_node, f_parts, False, s_type[1] if s_type is not None else None)

            if s_type is None and f_type is not None:
                s_type = self.get_tag_and_type_from_node(s_node, s_parts, True, f_type[1])

            if s_type is None and f_type is None:
                s_type, f_type = self.guess_connection_between_nodes(s_node, f_node)

            if s_type is None or f_type is None:
                raise Exception(f"Could not identify connection for edge {s} → {f}")

        # NOW: Let's build the connection
        s_node.connection_map[s_type[0]] = f_type[0], f_node

        correct_type = f_type[1].can_receive_from(s_type[1])

        if not correct_type:
            Logger.log(f"Mismatch of types when joining '{s}' to {f}' "
                       f"({s_type[1].id()} -/→ {f_type[1].id()})",
                       LogLevel.CRITICAL)
            Logger.log(f"No action taken to correct type-mismatch of '{s}' "
                       f"to {f}'")

        col = 'black' if correct_type else 'r'
        self.graph.add_edge(s_node, f_node, type_match=correct_type, color=col)

    @staticmethod
    def get_tag_and_type_from_node(node: Node, input_parts: List[str], is_start: bool, guess_type: DataType=None) -> Tuple[str, DataType]:
        lbl = input_parts[0]
        s = "/".join(input_parts)
        if node.node_type == NodeTypes.INPUT:
            if len(input_parts) != 1:
                message = f"The input tag '{s}' over-represents the INPUT node, this was corrected ({s} → {lbl}"
                Logger.log(message, LogLevel.WARNING)
            return lbl, next(iter(node.outputs().values())).output_type

        elif node.node_type == NodeTypes.OUTPUT:
            if len(input_parts) != 1:
                message = f"The input tag '{s}' over-represents the OUTPUT node, this was corrected ({s} → {lbl}"
                Logger.log(message, LogLevel.WARNING)
            return lbl, None

        elif node.node_type == NodeTypes.TASK:
            step_node: StepNode = node
            ins, outs = node.inputs(), node.outputs()
            put_name = 'output' if is_start else 'input'

            types = [(x, outs[x].output_type) for x in outs] \
                     if is_start else [(x, ins[x].input_type) for x in ins]
            if len(types) == 0:
                raise InvalidStepsException(f"The step '{input_parts}' referenced the tool "
                                            f"'{step_node.step.get_tool().id()}' with "
                                            f"no {put_name}s")

            if len(input_parts) == 1:
                # No tag, better hope there's only one input
                if len(types) > 1:
                    if guess_type is None:
                        Logger.log("Could not correctly determine the type for '{}', "
                                   "however it may be able to guess it later", LogLevel.WARNING)
                    else:
                        # try to get it here
                        raise NotImplementedError("The 'auto-align step by datatype' has not yet been implemented")

                else:
                    Logger.log(f"The tag '{s}' under-represents the STEP node, this was corrected "
                               f"({s} → {s}/{types[0][0]})", LogLevel.WARNING)

                    return types[0]
            else:
                # We have at least two parts, let's correct this
                tag = input_parts[1]
                t = [x for x in types if x[0] == tag]
                if len(t) == 0:
                    raise Exception(f"Could not identify an {put_name} called '{tag}' on the node '{node.id()}'")
                elif len(t) > 1:
                    raise Exception(f"Could not uniquely identify {put_name} with tag '{tag}'")
                else:
                    return t[0]

        raise Exception("Unhandled pathway in 'get_tag_and_type_from_node'")


    @staticmethod
    def guess_connection_between_nodes(s_node: Node, f_node: Node) -> Tuple[Tuple[str, DataType], Tuple[str, DataType]]:
        outs, ins = s_node.outputs(), f_node.inputs()

        s_types: List[Tuple[str, DataType]] = [(x, outs[x].output_type) for x in outs]
        f_types: List[Tuple[str, DataType]] = [(x, ins[x].input_type) for x in ins]

        # O(n**2) for determining types
        matching_types: List[Tuple[Tuple[str, DataType], Tuple[str, DataType]]] = []
        for s_type in s_types:
            matching_types.extend([(s_type, f_type) for f_type in f_types if f_type[1].can_receive_from(s_type)])

        if len(matching_types) == 0:
            raise InvalidInputsException(f"Can't find a common type between '{s_node.id()}' and '{f_node.id()}'")

        if len(matching_types) > 1:
            raise InvalidInputsException(f"Couldn't guess the single connection between '{s_node.id()}' "
                                         f"and '{f_node.id()}', there was more than 1 compatible connection")

        return matching_types[0]


    def attempt_connect_node(self, node: Node):
        """
        Get the inputs and outputs of the node, and try to make the connections
        :param node:
        """
        unfulfilled_connections: Dict[str, List[NodeAndTag]] = {}

        inputs: Dict[str, ToolInput] = node.inputs()
        outputs = node.outputs()

        if node.node_type == NodeTypes.INPUT:
            pass
            # try to fulfill unfulfilled connections

        if node.node_type == NodeTypes.TASK:
            pass










