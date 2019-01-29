from abc import ABC, abstractmethod
from typing import List

from janis import ToolOutput, ToolInput, String, File
from janis.tool.tool import Tool, ToolTypes
from janis.translations.cwl.cwl import Cwl


class ExpressionTool(Tool, ABC):

    def id(self) -> str:
        return self.tool()

    @staticmethod
    @abstractmethod
    def tool():
        raise NotImplementedError("override 'tool' method")

    @abstractmethod
    def expression(self):
        raise NotImplementedError("Must implement 'expression' method")

    def cwl(self):
        CLT = Cwl.ExpressionTool
        d = {
            Cwl.kCLASS: Cwl.Class.kCOMMANDLINETOOL,
            Cwl.kCWL_VERSION: Cwl.kCUR_VERSION,
            CLT.kID: self.id(),
            CLT.kLABEL: self.id(),
            CLT.kEXPRESSION: self.expression()
        }

        if self.inputs():
            inputs = self.inputs()
            for i in inputs:
                i.position = None
                i.prefix = None

            d[CLT.kINPUTS] = {t.tag: t.cwl() for t in inputs}

        if self.outputs():
            d[CLT.kOUTPUTS] = {t.tag: t.cwl() for t in self.outputs()}

        return d

    def type(self):
        return ToolTypes.ExpressionTool
