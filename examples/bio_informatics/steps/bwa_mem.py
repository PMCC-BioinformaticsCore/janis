from typing import List

from examples.bio_informatics.data_types.fastq import FastQ
from examples.bio_informatics.data_types.sam import Sam
from pipeline_definition.types.common_data_types import String, Number, Array
from pipeline_definition.types.tool import Tool, ToolOutput, ToolInput


class BwaMem(Tool):

    @staticmethod
    def tool():
        return "bwa-mem"

    @staticmethod
    def base_command():
        return "javac"

    def inputs(self) -> List[ToolInput]:
        return [
            ToolInput("reads", Array()),
            ToolInput("reference", FastQ()),
            ToolInput("outputFilename", String()),
            ToolInput("readGroup", String()),
            ToolInput("threads", Number(optional=True)),
            ToolInput("min_std_max_min", Array(optional=True))
        ]

    def outputs(self) -> List[ToolOutput]:
        return [
            ToolOutput("out", Sam())
        ]
