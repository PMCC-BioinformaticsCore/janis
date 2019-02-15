#
# Untar a file

from janis import Array, ToolInput, ToolOutput, InputSelector, File, Filename
from janis.unix.tools.unixtool import UnixTool


class Tar(UnixTool):

    @staticmethod
    def tool():
        return "Tar"

    def friendly_name(self):
        return "Tar archive"

    @staticmethod
    def base_command():
        return ["tar", "cvf"]

    def inputs(self):
        return [
            ToolInput("files", Array(File()), position=2),
            ToolInput("files2", Array(File(), optional=True), position=3),
            ToolInput("tarName", Filename(extension=".tar"), position=1)
        ]

    def outputs(self):
        return [
            ToolOutput("tarred", File(), glob=InputSelector("tarName"))     # "$(inputs.tarName)")
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