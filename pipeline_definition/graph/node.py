"""
    node.py

    Provides base class that different nodes must override, this translates closest to a Step
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple

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

    def __init__(self, node_type: int, label: str, depth=0):

        if label in self._N_node_map:
            raise Exception(f"Label {label} has already been used by node: {self._N_node_map[label]}")

        self.node_type: int = node_type
        self.label: str = label
        self.depth = depth

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

    def set_depth(self, depth: int):
        self.depth = max(self.depth, depth)

    @abstractmethod
    def inputs(self) -> Dict[str, ToolInput]:
        raise Exception(f"Subclass {type(self)} must implement inputs, return dict: key: StepInput")

    @abstractmethod
    def outputs(self) -> Dict[str, ToolOutput]:
        raise Exception(f"Subclass {type(self)} must implement outputs, return dict: key: StepOutput")


def layout_nodes(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    """
    Stack on depth away from root, and scale smaller columns
    :param nodes: list of nodes in Graph to position
    :return: Dict[Node, (x,y)]
    """
    pos = {}
    cur_depth = []
    depth_node: Dict[int, List[Node]] = {}

    def get_idx_for_depth(depth: int):
        d_idx = depth
        m = len(cur_depth)
        if d_idx >= m:
            cur_depth.extend([0] * (depth - m + 1))
        return cur_depth[depth]

    def push_to_dict(key, val):
        if key in depth_node:
            depth_node[key].append(val)
        else:
            depth_node[key] = [val]

    for node in nodes:
        depth = node.depth
        push_to_dict(depth, node)
        d = get_idx_for_depth(depth)
        cur_depth[depth] += 1
        pos[node] = (depth, d)

    # Now normalise each depth
    max_in_col = float(max(cur_depth))

    for (idx, d) in enumerate(cur_depth):
        if idx not in depth_node:
            continue
        nodes = depth_node[idx]
        scale = max_in_col / len(nodes)

        # If max and this col differ in even / odd, add half a unit
        half_bias = 0 if (len(nodes) % 2) == (max_in_col % 2) else 0.5

        for node in nodes:
            (x, y) = pos[node]
            pos[node] = (x, (y + half_bias)*scale)

    return pos


def layout_nodes2(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    # Aim for something like Rabix
    inputs = [n for n in nodes if n.node_type == NodeType.INPUT]
    others = [n for n in nodes if n.node_type != NodeType.INPUT]

    pos = layout_nodes(others)
    s = 0
    for n in inputs:
        pos[n] = (0, s)
        s += 1
    return pos


def layout_nodes3(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    # Other ideas to scale might be, start with steps and
    # stack in a similar way to depth around center, then do the same
    return None