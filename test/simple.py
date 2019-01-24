import yaml, os
from wehi.spec import Wehi
from Pipeline import Logger, Workflow

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
        translator = Wehi(label)
        translator.parse_string(_yml)

        # SimplePipeline.dump_cwl(translator.workflow, to_disk=False)
        SimplePipeline.dump_wdl(translator.workflow, to_disk=False)

        Logger.close_file()

    @staticmethod
    def dump_cwl(translator: Workflow, to_disk: False):
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
            with open(d + translator.id() + ".cwl", "w+") as cwl:
                Logger.log(f"Writing {translator.id()}.cwl to disk")
                yaml.dump(cwl_data, cwl, default_flow_style=False)
                Logger.log(f"Written {translator.id()}.cwl to disk")

            with open(d + translator.id() + "-job.yml", "w+") as cwl:
                Logger.log(f"Writing {translator.id()}-job.yml to disk")
                yaml.dump(inp_data, cwl, default_flow_style=False)
                Logger.log(f"Written {translator.id()}-job.yml to disk")

            for tool in tools_ar:
                tool_name = tool["id"].lower()
                with open(d_tools + tool_name + ".cwl", "w+") as cwl:
                    Logger.log(f"Writing {tool_name}.cwl to disk")
                    yaml.dump(tool, cwl, default_flow_style=False)
                    Logger.log(f"Written {tool_name}.cwl to disk")

    @staticmethod
    def dump_wdl(translator: Workflow, to_disk: False):
        import json
        wdl_data, inp_data, tools_dict = translator.wdl()
        print(wdl_data)
        print("================")
        print(inp_data)
        print("================")
        print("\n*******\n".join(tools_dict.values()))

        d = os.path.expanduser("~") + f"/Desktop/{translator.identifier}/wdl/"
        d_tools = d + "tools/"
        if not os.path.isdir(d):
            os.makedirs(d)
        if not os.path.isdir(d_tools):
            os.makedirs(d_tools)

        if to_disk:
            with open(d + translator.identifier + ".wdl", "w+") as wdl:
                Logger.log(f"Writing {translator.identifier}.wdl to disk")
                wdl.write(wdl_data)
                Logger.log(f"Written {translator.identifier}.wdl to disk")

            with open(d + translator.identifier + "-job.json", "w+") as inp:
                Logger.log(f"Writing {translator.identifier}-job.json to disk")
                json.dump(inp_data, inp)
                Logger.log(f"Written {translator.identifier}-job.json to disk")

            for tool_name in tools_dict:
                tool = tools_dict[tool_name]
                with open(d_tools + tool_name + ".wdl", "w+") as wdl:
                    Logger.log(f"Writing {tool_name}.cwl to disk")
                    wdl.write(tool)
                    Logger.log(f"Written {tool_name}.cwl to disk")


if __name__ == '__main__':
    SimplePipeline.test_simple()
