from abc import ABC, abstractmethod
import re
from typing import List, Dict, Optional, Any

from Pipeline.tool.tool import Tool, ToolArgument, ToolInput, ToolOutput, ToolTypes


class CommandTool(Tool, ABC):
    """
    Notes:
        - If you're thinking about secondary files, DON'T!! Consider creating a new DataType.
        - This class is similar to how RABIX COMPOSER creates the tools
        - You can subclass and override whichever fields you'd like, including the INPUTS / OUTPUTS!
        - Take note which options you can provide to the ToolInput and ToolOutput.
    """

    def id(self):
        return self.tool()

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

    def memory(self) -> Optional[int]:
        return None

    def cpus(self) -> Optional[int]:
        return None

    def arguments(self) -> Optional[List[ToolArgument]]:
        return None

    @classmethod
    def type(cls):
        return ToolTypes.CommandTool

    def inputs(self) -> List[ToolInput]:
        cls = type(self)
        attrs = [x for x in vars(cls) if not x.startswith("_")]
        d = cls.__dict__
        return [d[x] for x in attrs if x in d and issubclass(type(d[x]), ToolInput)]

    def outputs(self) -> List[ToolOutput]:
        cls = type(self)
        attrs = [x for x in vars(cls) if not x.startswith("_")]
        d = cls.__dict__
        return [d[x] for x in attrs if x in attrs and issubclass(type(d[x]), ToolOutput)]

    def cwl(self) -> Dict[str, Any]:
        d = {
            "class": "CommandLineTool",
            "cwlVersion": "v1.0",
            "baseCommand": self.base_command(),
            "id": self.tool(),
            "label": self.tool(),
        }

        hints = {}
        if self.docker() is not None:
            hints["DockerRequirement"] = {"dockerPull": self.docker()}

        if hints:
            d["hints"] = hints

        inps = {}
        for tool_input in self.inputs():
            inps[tool_input.tag] = tool_input.cwl()

        if self.inputs():
            d["inputs"] = {t.tag: t.cwl() for t in self.inputs()}

        if self.outputs():
            d["outputs"] = {t.tag: t.cwl() for t in self.outputs()}

        args = self.arguments()
        if args and args is not None:
            d["arguments"] = [a.cwl() for a in args]

        return d

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
