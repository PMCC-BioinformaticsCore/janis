import examples.unix_commands
from pipeline_definition.utils.logger import Logger, LogLevel
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


class SimplePipeline:

    @staticmethod
    def test_simple():
        label = "simple"
        Logger.set_write_location("/Users/franklinmichael/source/wehi-pipeline-definition/simple.log")
        translator = PipelineTranslator(label)
        translator.translate_string(_yml)

        # SimplePipeline.dump_cwl(translator)
        SimplePipeline.dump_wdl(translator)

        Logger.close_file()

    @staticmethod
    def dump_cwl(translator: PipelineTranslator):
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

        with open(d + translator.label + ".cwl", "w+") as cwl:
            Logger.log(f"Writing {label}.cwl to disk")
            yaml.dump(cwl_data, cwl, default_flow_style=False)
            Logger.log(f"Written {label}.cwl to disk")

        with open(d + translator.label + "-job.yml", "w+") as cwl:
            Logger.log(f"Writing {label}-job.yml to disk")
            yaml.dump(inp_data, cwl, default_flow_style=False)
            Logger.log(f"Written {label}-job.yml to disk")

        for tool in tools_ar:
            tool_name = tool["id"].lower()
            with open(d_tools + tool_name + ".cwl", "w+") as cwl:
                Logger.log(f"Writing {tool_name}.cwl to disk")
                yaml.dump(tool, cwl, default_flow_style=False)
                Logger.log(f"Written {tool_name}.cwl to disk")

    @staticmethod
    def dump_wdl(translator: PipelineTranslator):
        wdl_data, inp_data, tools_ar = translator.wdl()
        print(wdl_data)
        print("================")
        print(inp_data)
        print("================")
        print("*******".join(tools_ar))

if __name__ == '__main__':
    SimplePipeline.test_simple()
