from typing import List

from pipeline_definition.types.tool import Tool, ToolOutput, ToolInput


class BwaMem(Tool):
    @staticmethod
    def tool():
        return "bwa-mem"

    @staticmethod
    def supported_translations() -> List[str]:
        return ["cwl"]

    def inputs(self) -> List[ToolInput]:
        return [ToolInput()]

    def outputs(self) -> List[ToolOutput]:
        pass