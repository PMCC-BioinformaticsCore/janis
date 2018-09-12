# https://github.com/common-workflow-language/common-workflow-language/blob/master/v1.0/examples/arguments.cwl

cwlVersion: v1.0
class: CommandLineTool
label: Example trivial wrapper for Java 7 compiler
hints:
  DockerRequirement:
    dockerPull: java:7-jdk
baseCommand: javac
arguments: ["-d", $(runtime.outdir)]
inputs:
  src:
    type: File
    inputBinding:
      position: 1
outputs:
  classfile:
    type: File
    outputBinding:
      glob: "*.class"