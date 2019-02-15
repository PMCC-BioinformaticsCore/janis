"""
    node.py

    Provides base class that different nodes must override, this translates closest to a Step
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Any

from janis.tool.tool import ToolInput, ToolOutput

NodeLabel = str
NodeType = int


class NodeTypes:
    INPUT: NodeType = 1
    OUTPUT: NodeType = 2
    TASK: NodeType = 3

    @staticmethod
    def to_str(node_type: NodeType) -> str:
        if node_type == NodeTypes.INPUT: return "Input"
        if node_type == NodeTypes.OUTPUT: return "Output"
        if node_type == NodeTypes.TASK: return "Task"
        raise Exception(f"Unhandled task type: '{node_type}'")

    @staticmethod
    def to_col(node_type: NodeType) -> str:
        if node_type == NodeTypes.INPUT: return "red"
        if node_type == NodeTypes.OUTPUT: return "lightblue"
        if node_type == NodeTypes.TASK: return "blue"
        raise Exception(f"Unhandled task type: '{node_type}'")


class Node(ABC):

    _N_counter: int = 1
    _N_nodeId_map: Dict[int, Any] = {}

    def __init__(self, node_type: NodeType, label: NodeLabel, depth=0):

        self.node_type: NodeType = node_type
        self._label: NodeLabel = label
        self.depth = depth

        self.connection_map: Dict[str, Any] = {}    # actually an edge

        # Update unique counter for hash
        self._nodeId = Node._N_counter
        Node._N_counter += 1

        # Map the node, so we can look it up later
        self._N_nodeId_map[self._nodeId] = self

    def id(self) -> str:
        return self._label

    def __hash__(self):
        return self._nodeId

    def __repr__(self):
        return f"{self.node_type}: {self.id()}"

    def __str__(self):
        return f"{NodeTypes.to_str(self.node_type)}: {self._label}"

    def set_depth(self, depth: int):
        self.depth = max(self.depth, depth)

    def __eq__(self, other):
        if isinstance(other, Node):
            return self._nodeId == other._nodeId
        return False

    @abstractmethod
    def inputs(self) -> Dict[str, ToolInput]:
        raise Exception(f"Subclass {type(self)} must implement inputs, return dict: key: ToolInput")

    @abstractmethod
    def outputs(self) -> Dict[str, ToolOutput]:
        raise Exception(f"Subclass {type(self)} must implement outputs, return dict: key: ToolOutput")

    def __setitem__(self, key, value):
        self.connection_map[key] = value


NodeAndTag = Tuple[str, Node]


def layout_nodes(nodes: List[Node], n_inputs: int=0) -> Dict[Node, Tuple[int, int]]:
    """
    Stack on depth away from root, and scale smaller columns
    :param n_inputs:
    :param nodes: list of nodes in Graph to position
    :return: Dict[Node, (x,y)]
    """
    pos = {}
    # int is the depth at the specified index
    cur_depth: List[int] = []
    depth_node: Dict[int, List[Node]] = {}

    def get_idx_for_depth(d: int):
        d_idx = d
        m = len(cur_depth)
        if d_idx >= m:
            cur_depth.extend([0] * (d - m + 1))
        return cur_depth[d]

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
    max_in_col = float(max(cur_depth + [n_inputs]))

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
    inputs = [n for n in nodes if n.node_type == NodeTypes.INPUT]
    others = [n for n in nodes if n.node_type != NodeTypes.INPUT]

    pos = layout_nodes(others, len(inputs))
    s = 0
    for n in inputs:
        pos[n] = (0, s)
        s += 1
    return pos


def layout_nodes3(nodes: List[Node]) -> Dict[Node, Tuple[int, int]]:
    # Other ideas to scale might be, start with steps and
    # stack in a similar way to depth around center, then do the same
    return {}