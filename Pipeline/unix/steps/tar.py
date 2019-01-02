#
# Untar a file
from typing import List

from Pipeline import Array
from Pipeline.types.common_data_types import File, Filename
from Pipeline.tool.commandtool import CommandTool, ToolInput, ToolOutput


class Tar(CommandTool):

    @staticmethod
    def tool():
        return "Tar"

    @staticmethod
    def base_command():
        return ["tar", "cvf"]

    def inputs(self):
        return [
            ToolInput("files", Array(File()), position=2),
            ToolInput("tarName", Filename(extension=".tar"), position=1)
        ]

    def outputs(self):
        return [
            ToolOutput("tarred", File(), glob="*.tar")
        ]


if __name__ == "__main__":
    print(Tar().help())

"""

EXAMPLE XML

        <Step>
            <title>Tar</title>
            <author>Michael Franklin</author>
            <cores>1</cores>
            <ram>1024</ram>
            <label>{{label}}</label>
            <supportedTypes>
                <type>cwl</type>
                <type>wdl</type>
            </supportedTypes>
        
            <inputs>
                <input>
                    <tag>input2</tag>
                    <type>File</type>
                </input>
                <input>
                    <tag>input2</tag>
                    <type>File</type>
                </input>
            </inputs>
            <outputs>
                <output>
                    <tag>output</tag>
                    <type>File</type>
                </output>
            </outputs>
        </Step>

==============================

        class: CommandLineTool
        cwlVersion: v1.0
        id: tar
        baseCommand:
          - tar
          - cvf
        inputs:
          - id: tarName
            type: string?
            inputBinding:
              position: 0
          - id: input1
            type:
              - File
            inputBinding:
              position: 1
          - id: input2
            type:
              - File
            inputBinding:
              position: 2
        outputs:
          - id: output
            type: File?
            outputBinding:
              glob: '*'
        label: tar

  ========================================
  
        task SingleTar {
            String tarName
            File input1
            File input2
        
            command {
                tar cvf ${tarName} ${input1} ${input2}
            }
        
            output {
                File out = glob("*")[0]
            }
        }
"""