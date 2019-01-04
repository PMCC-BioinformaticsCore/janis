from abc import ABC, abstractmethod
import re
from typing import List, Dict, Optional, Any

from Pipeline.tool.tool import Tool, ToolArgument, ToolInput, ToolOutput, ToolTypes

import cwlgen.cwlgen as cwl


class CommandTool(Tool, ABC):
    """
    Notes:
        - If you're thinking about secondary files, DON'T!! Consider creating a new DataType.
        - This class is similar to how RABIX COMPOSER creates the tools
        - You can subclass and override whichever fields you'd like, including the INPUTS / OUTPUTS!
        - Take note which options you can provide to the ToolInput and ToolOutput.
    """

    def id(self):
        return self.tool().lower()

    @staticmethod
    @abstractmethod
    def tool():
        raise Exception(f"subclass MUST implement 'tool' method")

    @classmethod
    def full_name(cls):
        if cls.version() is not None:
            return f"{cls.tool()}/{cls.version()}"
        return cls.tool()

    @staticmethod
    def version():
        return None

    @staticmethod
    def docker():
        return "ubuntu:latest"

    @staticmethod
    @abstractmethod
    def base_command():
        raise Exception("Subclass MUST implement the 'base_command' method with str or [str] return types")

    @staticmethod
    def stdout() -> Optional[str]:
        return None

    def memory(self) -> Optional[int]:
        return None

    def cpus(self) -> Optional[int]:
        return None

    def arguments(self) -> Optional[List[ToolArgument]]:
        return None

    @classmethod
    def type(cls):
        return ToolTypes.CommandTool

    @staticmethod
    def environment_variables() -> Optional[Dict[str, str]]:
        return None

    def cwl(self, with_docker=True) -> Dict[str, Any]:
        tool = cwl.CommandLineTool(
            tool_id=self.id(),
            base_command=self.base_command(),
            label=self.id(),
            doc=self.doc(),
            # cwl_version=Cwl.kCUR_VERSION,
            stdin=None,
            stderr=None,
            stdout=self.stdout()
        )

        tool.requirements.extend([
            cwl.InlineJavascriptReq()
        ])

        if with_docker:
            tool.requirements.append(cwl.DockerRequirement(
                docker_pull=self.docker(),
                # docker_load=None,
                # docker_file=None,
                # docker_import=None,
                # docker_image_id=None,
                # docker_output_dir=None
            ))

        tool.inputs.extend(i.cwl() for i in self.inputs())
        tool.outputs.extend(o.cwl() for o in self.outputs())
        args = self.arguments()
        if args is not None:
            tool.arguments.extend(a.cwl() for a in self.arguments())

        return tool.get_dict()

    def _command(self):

        base = self.base_command() if type(self.base_command()) == str else " ".join(self.base_command())
        command = base

        if self.arguments():
            args: List[ToolArgument] = sorted(self.arguments(), key=lambda a: a.position)
            args_joined = " ".join(a.wdl() for a in args)
            base += " " + args_joined

        if self.inputs():
            inps = " ".join(f"${{{s.tag}}}" for s in sorted(self.inputs(), key=lambda a: a.position))
            command += " " + inps

        return command

    def wdl_name(self):
        return self.tool().replace("-", "_")

    def wdl(self):

        command = self._command()
        tb = "    "
        nl = "\n"

        inputs_str = "{tb}{data_type} {identifier}{default_with_equals_if_required}"
        outputs_str = "{tb2}{data_type} {identifier} = glob(\"{glob}\")[0]"

        default_with_equals = lambda d: f' = "{str(d)}"' if d is not None else ""
        inputs = nl.join([inputs_str.format(
            tb=tb,
            data_type=i.input_type.wdl(),
            identifier=i.tag,
            default_with_equals_if_required=default_with_equals(i.input_type.default())
        ) for i in self.inputs()])
        outputs = nl.join([outputs_str.format(
            tb2=2 * tb,
            data_type=o.output_type.wdl(),
            identifier=o.tag,
            glob=o.glob
        ) for o in self.outputs()])

        # inputs = "\n\t".join([f"{t.input_type.wdl()} {t.tag}" for t in self.inputs()])

        # outputs = "\n\t\t".join(f"{o.output_type.wdl()} {o.tag} = glob(\"{o.glob}\")[0]" for o in self.outputs())

        return f"""
task {self.wdl_name()} {{
{inputs}
    
    runtime {{ docker: "{self.docker()}" }}
    command {{
        {command}
    }}
    
    output {{
{outputs}
    }}
}}"""

    def help(self):
        import inspect
        path = inspect.getfile(self.__class__)

        ins = sorted(self.inputs(), key=lambda i: i.position if i.position is not None else 0)
        args = ""
        if self.arguments():
            args = " " + " ".join(f"{(a.prefix if a.prefix is not None else '') + ' ' if (a.prefix is not None and a.separate_value_from_prefix) else ''}{a.value}" for a in self.arguments())

        prefixes = " -" + "".join(i.prefix.replace("-", "").replace(" ", "") for i in ins if i.prefix is not None)

        docker = self.docker()

        command = (self.base_command() if isinstance(self.base_command(), str) else " ".join(self.base_command())) \
                  + args + prefixes

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (t.prefix + ": ") if t.separate_value_from_prefix else t.prefix
            return f"\t\t{t.tag} ({prefix_with_space}{t.input_type.id()}{('=' + str(t.default)) if t.default is not None else ''})" \
                f": {'' if t.doc is None else t.doc}"

        output_format = lambda t: f"\t\t{t.tag} ({t.output_type.id()}): {'' if t.doc is None else t.doc}"

        requiredInputs = "\n".join(input_format(x) for x in ins if not x.optional)
        optionalInputs = "\n".join(input_format(x) for x in ins if x.optional)
        outputs = "\n".join(output_format(o) for o in self.outputs())

        return f"""
    Pipeline tool: {path} ({self.id()})
NAME
    {self.id()}

SYNOPSIS
    {command}

DOCKER
    {docker}

DESCRIPTION
    {self.doc() if self.doc is not None else "No documentation provided"}

INPUTS:
    REQUIRED:
{requiredInputs}

    OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
"""
