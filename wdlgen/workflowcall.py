from abc import ABC, abstractmethod
from typing import Dict, List

from wdlgen.task import Task
from wdlgen.util import WdlBase


class WorkflowCallBase(WdlBase, ABC):
    @abstractmethod
    def get_string(self, indent=1):
        raise Exception("Must override 'get_string(indent:int)'")


class WorkflowCall(WorkflowCallBase):

    def __init__(self, namespaced_identifier: str=None, alias: str=None, inputs_map: Dict[str, str]=None):
        """

        :param task:
        :param namespaced_identifier: Required if task is imported. The workflow might take care of this later?
        :param alias:
        :param inputs_map:
        """
        self.namespaced_identifier = namespaced_identifier
        self.alias = alias
        self.inputs_map = inputs_map

        self.call_format = """{ind}call {name}{alias} {body}"""
        self.body_format = """{{\n{ind}{tb}input:\n{input_map}\n{ind}}}"""

    def get_string(self, indent=1):
        tb = "  "

        body = ""
        if self.inputs_map:
            inpmap = ",\n".join(((indent + 2) * tb + str(k) + "=" + str(v)) for k, v in self.inputs_map.items())
            body = self.body_format.format(ind=indent * tb, input_map=inpmap, tb=tb)

        return self.call_format.format(
            ind=indent * tb,
            tb=tb,
            name=self.namespaced_identifier if self.namespaced_identifier else self.task.name,
            alias=(" as " + self.alias) if self.alias else "",
            body=body
        )


class WorkflowScatter(WorkflowCallBase):

    def __init__(self, identifier: str, expression: str, calls: List[WorkflowCall]=None):
        self.identifier: str = identifier
        self.expression: str = expression
        self.calls: List[WorkflowCall] = calls if calls else []

    def get_string(self, indent=1):
        scatter_iteration_statement = "{identifier} in {expression}"\
            .format(identifier=self.identifier, expression=self.expression)

        body = "\n".join(c.get_string(indent=indent + 1) for c in self.calls)

        return "{ind}scatter ({st}) {{\n {body}\n{ind}}}"\
            .format(ind=indent * '  ', st=scatter_iteration_statement, body=body)
