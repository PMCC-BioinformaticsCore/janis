from typing import Dict, List, Tuple, Set

from Pipeline.types.data_types import DataType
from Pipeline.utils.logger import Logger, LogLevel
from Pipeline.utils.descriptor import BaseDescriptor, GetOnlyDescriptor

import networkx as nx

from Pipeline.graph.node import Node, NodeTypes, NodeTag
from Pipeline.utils.errors import DuplicateLabelIdentifier, InvalidNodeIdentifier, NodeNotFound
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
        self.connections: Dict[str, Set[NodeTag]] = {}

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

        if s_node.node_type == NodeTypes.INPUT:
            if len(s_parts) != 1:
                message = f"The input tag '{s}' over-represents the INPUT node, this was corrected ({s} → {s_label}"
                Logger.log(message, LogLevel.WARNING)
        elif s_node.node_type ==


    @staticmethod
    def get_type_from_node(node: Node, input_parts: List[str], is_start: bool) -> Tuple[str, DataType]:
        lbl = input_parts[0]
        if node.node_type == NodeTypes.INPUT:
            if len(input_parts) != 1:
                message = f"The input tag '{s}' over-represents the INPUT node, this was corrected ({s} → {lbl}"
                Logger.log(message, LogLevel.WARNING)
            return lbl, next(iter(node.outputs().values())).output_type

        elif node.node_type == NodeTypes.OUTPUT:
            if len(input_parts) != 1:
                message = f"The input tag '{s}' over-represents the OUTPUT node, this was corrected ({s} → {lbl}"
                Logger.log(message, LogLevel.WARNING)
            return lbl, next(iter(node.inputs().values())).input_type

        elif node.node_type == NodeTypes.TASK:
            ins, outs = node.inputs(), node.outputs()
            types = [outs[x].output_type for x in outs] \
                     if is_start else [ins[x].input_type for x in ins]
            if len(input_parts) == 1:
                # No tag, better hope there's only one input



    def attempt_connect_node(self, node: Node):
        """
        Get the inputs and outputs of the node, and try to make the connections
        :param node:
        """
        unfulfilled_connections: Dict[str, List[NodeTag]] = {}

        inputs: Dict[str, ToolInput] = node.inputs()
        outputs = node.outputs()

        if node.node_type == NodeTypes.INPUT:
            pass
            # try to fulfill unfulfilled connections

        if node.node_type == NodeTypes.TASK:
            pass










