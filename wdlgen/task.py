from typing import List, Optional

from wdlgen.types import WdlType
from wdlgen.util import WdlBase


class Runtime(WdlBase):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def wdl(self):
        return ["{k}: {v}".format(k=k, v=v) for k,v in self.kwargs.items()]


class Input(WdlBase):
    def __init__(self, data_type: WdlType, name: str, expression: str=None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        return "{type} {name}{def_w_equals}".format(
            type=self.type.wdl(),
            name=self.name,
            def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
        )


class Output(WdlBase):
    def __init__(self, data_type: WdlType, name: str, expression: str=None):
        self.type = data_type
        self.name = name
        self.expression = expression

    def wdl(self):
        return "{type} {name}{def_w_equals}".format(
            type=self.type.wdl(),
            name=self.name,
            def_w_equals=(" = {val}".format(val=self.expression) if self.expression else "")
        )


class Command(WdlBase):
    """
    Past the regular attributes, I've built the command generation here, because that's where
    it logically is. This logic is pretty similar to CWL (and Rabix Composer's visualisation).
    Declare a base command, arguments, inputs with positions and prefixes and we can manually assemble the command.

    For this reason, there's a bit of double handling, if you add an input to a Task,
    you also need to add it here, and that's just unfortunate.

    ======

    Just want a plain and simple command, no worries just don't add inputs or anything else.

    """

    class CommandArgument:
        def __init__(self, prefix: str=None, position: int=None):
            self.prefix: str = prefix
            self.position: int = position

        def wdl(self):
            return self.prefix if self.prefix else ""

    class CommandInput(CommandArgument):
        def __init__(self, name: str, optional: bool=False, prefix: str=None, position: int=None, separate_value_from_prefix: bool=True):
            super().__init__(prefix, position)
            self.name = name
            self.optional = optional
            self.separate = separate_value_from_prefix

        @staticmethod
        def from_input(inp: Input, prefix: str=None, position: int=None):
            return Command.CommandInput(inp.name, inp.type._optional, prefix, position)

        def wdl(self):
            sp = " " if self.separate else ""
            pr = self.prefix if self.prefix else ""
            bc = pr + sp

            if self.optional:
                return '${{"{pre}" + {val}}}'.format(pre=bc, val=self.name)
            else:
                return bc + "${{{val}}}".format(val=self.name)

    def __init__(self, command, inputs: Optional[List[CommandInput]]=None, arguments: Optional[List[CommandArgument]]=None):
        self.command = command
        self.inputs = inputs if inputs else []
        self.arguments = arguments if arguments else []

    def wdl(self, indent:int=0):
        if not self.command:
            raise Exception("No base 'command' has been set on this command object")
        if not (self.inputs or self.arguments):
            return self.command

        # build up command
        args = sorted([*self.inputs, *self.arguments], key=lambda a: a.position if a.position else 0)
        command: str = self.command if isinstance(self.command, str) else " ".join(self.command)
        tb = "  "
        tbed_arg_indent = tb * (indent + 1)

        return indent * tb + command + "".join([" \\\n" + tbed_arg_indent + a.wdl() for a in args])


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

    def __init__(self, identifier: str, inputs: List[Input]=None, outputs: List[Output]=None, command: Command=None, runtime: Runtime=None):
        self.identifier = identifier
        self.inputs = inputs if inputs else []
        self.outputs = outputs if outputs else []
        self.command = command
        self.runtime = runtime

        self.format = """
task {identifier} {{
{inputs_block}
{command_block}
{runtime_block}
{output_block}
}}
        """.strip()

    def wdl(self):
        tb = "  "

        identifier = self.identifier
        inputs_block, command_block, runtime_block, output_block = "", "", "", ""

        if self.inputs:
            inputs_block = "\n".join(tb + i.wdl() for i in self.inputs)

        if self.outputs:
            output_block = "{tb}output {{\n{outs}\n{tb}}}".format(
                tb=tb,
                outs="\n".join((2*tb) + o.wdl() for o in self.outputs)
            )

        if self.command:
            command_block = "{tb}command {{\n{args}\n{tb}}}".format(
                tb=tb,
                args=self.command.wdl(indent=2)
            )

        if self.runtime:
            runtime_block = "{tb}runtime {{\n{args}\n{tb}}}".format(
                tb=tb,
                args="\n".join((2 * tb) + a for a in self.runtime.wdl())
            )

        return self.format.format(
            identifier=identifier,
            inputs_block=inputs_block,
            command_block=command_block,
            runtime_block=runtime_block,
            output_block=output_block
        )
