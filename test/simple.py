from utils.logger import Logger
from pipeline_definition.pipeline_translator import PipelineTranslator
import yaml, os

_yml = """
inputs:
    tarfile:
        type: TarFile
        path: /Users/franklinmichael/source/simple-workflow/hello.tar
    tarName: hello.tar

steps:
    untar:
        tool: untar
        input: tarfile

    compile:
        tool: java-compiler
        input: untar/outp

    tar:
        tool: tar
        tarName: tarName
        input1: untar/outp
        input2: compile/outp
        
outputs:
    untarred: untar/outp
    compiled: compile
    tarred: tar/outp
"""

wdl = """
import tools/untar.wdl as U
import tools/java-compiler.wdl as J
import tools/tar.wdl as T

workflow simple {
	File tarfile
	string tarName
	call U.Untar { <wdl-mapping here> }
	call J.java-compiler { <wdl-mapping here> }
	call T.Tar { <wdl-mapping here> }

    output {
	File untarred = u.untarred
	File compiled = c.compiled
	File tarred = t.tarred
    }
}
"""

cwl = """
class: Workflow
cwlVersion: v1.0
id: simple
inputs:
  tarName:
    type: string
  tarfile:
    type: File
label: simple
outputs:
  compiled:
    outputSource: compile/outp
    type: File
  tarred:
    outputSource: tar/outp
    type: File
  untarred:
    outputSource: untar/outp
    type: File?
requirements:
- class: InlineJavascriptRequirement
steps:
  compile:
    in:
      input: untar/outp
    label: compile
    out:
    - outp
    run: tools/java-compiler.cwl
  tar:
    in:
      input1: untar/outp
      input2: compile/outp
      tarName: tarName
    label: tar
    out:
    - outp
    run: tools/tar.cwl
  untar:
    in:
      input: tarfile
    label: untar
    out:
    - outp
    run: tools/untar.cwl
"""


class SimplePipeline:

    @staticmethod
    def test_simple():

        label = "simple"
        Logger.set_write_location("/Users/franklinmichael/source/wehi-pipeline-definition/simple.log")
        translator = PipelineTranslator(label)
        translator.translate_string(_yml)

        # SimplePipeline.dump_cwl(translator, to_disk=False)
        SimplePipeline.dump_wdl(translator, to_disk=False)

        Logger.close_file()

    @staticmethod
    def dump_cwl(translator: PipelineTranslator, to_disk: False):
        cwl_data, inp_data, tools_ar = translator.cwl()

        d = os.path.expanduser("~") + "/Desktop/cwl-test/"
        d_tools = d + "tools/"

        if not os.path.isdir(d):
            os.mkdir(d)
        if not os.path.isdir(d_tools):
            os.mkdir(d_tools)

        print(yaml.dump(cwl_data, default_flow_style=False))
        print(yaml.dump(inp_data, default_flow_style=False))
        print(yaml.dump(tools_ar, default_flow_style=False))

        if to_disk:
            with open(d + translator.label + ".cwl", "w+") as cwl:
                Logger.log(f"Writing {translator.label}.cwl to disk")
                yaml.dump(cwl_data, cwl, default_flow_style=False)
                Logger.log(f"Written {translator.label}.cwl to disk")

            with open(d + translator.label + "-job.yml", "w+") as cwl:
                Logger.log(f"Writing {translator.label}-job.yml to disk")
                yaml.dump(inp_data, cwl, default_flow_style=False)
                Logger.log(f"Written {translator.label}-job.yml to disk")

            for tool in tools_ar:
                tool_name = tool["id"].lower()
                with open(d_tools + tool_name + ".cwl", "w+") as cwl:
                    Logger.log(f"Writing {tool_name}.cwl to disk")
                    yaml.dump(tool, cwl, default_flow_style=False)
                    Logger.log(f"Written {tool_name}.cwl to disk")

    @staticmethod
    def dump_wdl(translator: PipelineTranslator, to_disk: False):
        wdl_data, inp_data, tools_ar = translator.wdl()
        print(wdl_data)
        print("================")
        print(inp_data)
        print("================")
        print("*******\n".join(tools_ar))


if __name__ == '__main__':
    SimplePipeline.test_simple()
