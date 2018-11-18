import examples.unix_commands
from pipeline_definition.utils.logger import Logger, LogLevel

_yml = """
inputs:
    tarfile:
        type: File
        path: /Users/franklinmichael/source/simple-workflow/hello.tar
    tarName: 
        type: File
        path: hello.tar

outputs:
    untarred: untar/out
    compiled: compile/out
    tarred: tar/out

steps:
    untar:
        tool: untar
        input: tarfile

    compile:
        tool: java-compiler
        input: untar/out

    tar:
        tool: tar
        tarName: tarName
        input1: untar/out
        input2: compile/out
"""

class SimplePipeline():

    @staticmethod
    def test_simple():
        from pipeline_definition.pipeline_translator import PipelineTranslator
        Logger.CONSOLE_LEVEL = LogLevel.DEBUG
        translator = PipelineTranslator()
        translation = translator.translate_string(_yml)
        print(translation)


if __name__ == '__main__':
    SimplePipeline.test_simple()
