"""
    node.py

    Provides base class that different nodes must override, this translates closest to a Step
"""

from typing import Optional, List

class NodeType:
    INPUT = 1
    OUTPUT = 2
    TASK = 3

    @staticmethod
    def to_str(node_type: int):
        if node_type == NodeType.INPUT: return "Input"
        if node_type == NodeType.OUTPUT: return "Output"
        if node_type == NodeType.TASK: return "Task"


class Node:
    def __init(self, node_type: int, label: str):
        self.node_type: int = node_type
        self.label: str = label

    def __hash__(self):
        return self.label

    def __str__(self):
        return f"NODE {NodeType.to_str(self.node_type)}: {self.label}"

    def inputs(self) -> None:
        pass

    def outputs(self) -> List[None]:
        pass
