from wdlgen.util import WdlBase


class Task(WdlBase):

    """
    A task definition is a way of encapsulating a UNIX command and environment
    and presenting them as functions. Tasks have both inputs and outputs.
    Inputs are declared as declarations at the top of the task definition,
    while outputs are defined in the output section.

    The user must provide a value for these two parameters in order for this task to be runnable.

    A task is a declarative construct with a focus on constructing a command from a template.
    The command specification is interpreted in an engine specific way, though a typical case
    is that a command is a UNIX command line which would be run in a Docker image.

    Tasks also define their outputs, which is essential for building dependencies between tasks.
    Any other data specified in the task definition (e.g. runtime information and meta-data) is optional.

    Documentation: https://github.com/openwdl/wdl/blob/master/versions/draft-2/SPEC.md#task-definition
    """

    def __init__(self, identifier: str, inputs, command, runtime, output):
        self.identifier = identifier
        self.inputs = inputs
        self.command = command
        self.runtime = runtime
        self.output = output

        self.format = """
task {identifier} {
{inputs_block}
{command_block}
{runtime_block}
{output_block}
}
        """.strip()

    def wdl(self):
        identifier = self.identifier
        inputs_block, command_block, runtime_block, output_block = "", "", "", ""

        self.format.format(
            identifier=identifier,
            inputs_block=inputs_block,
            command_block=command_block,
            runtime_block=runtime_block,
            output_block=output_block
        )


class Command(WdlBase):
    def __init__(self, command):
        self.command = command

    def wdl(self):
        return None


class Runtime(WdlBase):
    pass

class Input(WdlBase):
    pass

class Output(WdlBase):
    pass