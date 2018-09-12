# https://github.com/common-workflow-language/common-workflow-language/blob/master/v1.0/examples/tar-param.cwl

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [tar, xf]
inputs:
  tarfile:
    type: File
    inputBinding:
      position: 1
  extractfile:
    type: string
    inputBinding:
      position: 2
outputs:
  example_out:
    type: File
    outputBinding:
      glob: $(inputs.extractfile)
