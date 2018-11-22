import examples.unix_commands
from pipeline_definition.utils.logger import Logger, LogLevel

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
        import yaml, os
        from pipeline_definition.pipeline_translator import PipelineTranslator
        label = "simple"
        Logger.set_write_location("/Users/franklinmichael/source/wehi-pipeline-definition/simple.log")
        translator = PipelineTranslator(label)
        translator.translate_string(_yml)
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

        with open(d + label + ".cwl", "w+") as cwl:
            Logger.log(f"Writing {label}.cwl to disk")
            yaml.dump(cwl_data, cwl, default_flow_style=False)
            Logger.log(f"Written {label}.cwl to disk")

        with open(d + label + "-job.yml", "w+") as cwl:
            Logger.log(f"Writing {label}-job.yml to disk")
            yaml.dump(inp_data, cwl, default_flow_style=False)
            Logger.log(f"Written {label}-job.yml to disk")

        for tool in tools_ar:
            tool_name = tool["id"].lower()
            with open(d_tools + tool_name + ".cwl", "w+") as cwl:
                Logger.log(f"Writing {tool_name}.cwl to disk")
                yaml.dump(tool, cwl, default_flow_style=False)
                Logger.log(f"Written {tool_name}.cwl to disk")

        Logger.close_file()

if __name__ == '__main__':
    SimplePipeline.test_simple()
