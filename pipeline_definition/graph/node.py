"""
    node.py

    Provides base class that different nodes must override, this translates closest to a Step
"""
from abc import ABC, abstractmethod
from typing import Dict

from pipeline_definition.types.tool import ToolInput, ToolOutput


class NodeType:
    INPUT = 1
    OUTPUT = 2
    TASK = 3

    @staticmethod
    def to_str(node_type: int):
        if node_type == NodeType.INPUT: return "Input"
        if node_type == NodeType.OUTPUT: return "Output"
        if node_type == NodeType.TASK: return "Task"

    @staticmethod
    def to_col(node_type: int):
        if node_type == NodeType.INPUT: return "red"
        if node_type == NodeType.OUTPUT: return "lightblue"
        if node_type == NodeType.TASK: return "blue"


class Node(ABC):

    _N_counter = 1
    _N_nodeId_map = {}
    _N_node_map = {}

    def __init__(self, node_type: int, label: str):

        if label in self._N_node_map:
            raise Exception(f"Label {label} has already been used by node: {self._N_node_map[label]}")

        self.node_type: int = node_type
        self.label: str = label

        # Update unique counter for hash
        self._nodeId = Node._N_counter
        Node._N_counter += 1

        # Map the node, so we can look it up later
        self._N_nodeId_map[self._nodeId] = self
        self._N_node_map[self.label] = self

    def __hash__(self):
        return self._nodeId

    def __str__(self):
        return f"{NodeType.to_str(self.node_type)}: {self.label}"

    @abstractmethod
    def inputs(self) -> Dict[str, ToolInput]:
        raise Exception(f"Subclass {type(self)} must implement inputs, return dict: key: StepInput")

    @abstractmethod
    def outputs(self) -> Dict[str, ToolOutput]:
        raise Exception(f"Subclass {type(self)} must implement outputs, return dict: key: StepOutput")
