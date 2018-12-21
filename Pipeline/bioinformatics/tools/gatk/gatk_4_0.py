from abc import ABC

from Pipeline import CommandTool, ToolInput, Boolean
from Pipeline.bioinformatics.tools.gatk.gatkbase import GatkBase


class Gatk_4_0(GatkBase, ABC):

    @staticmethod
    def docker():
        return "4.0.12.0"
