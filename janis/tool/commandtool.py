from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any

import cwlgen as cwl

from janis.tool.tool import Tool, ToolArgument, ToolInput, ToolTypes
from janis.utils.metadata import Metadata


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

    @staticmethod
    def environment_variables() -> Optional[Dict[str, str]]:
        return None

    @staticmethod
    def requirements() -> Optional[List[cwl.Requirement]]:
        return None

    @staticmethod
    def hint_map() -> Optional[Dict[str, Any]]:
        return None

    def translate(self, translation, with_docker=True):
        import janis.translations
        return janis.translations.translate_tool(self, translation, with_docker=with_docker)

    def help(self):
        import inspect
        path = inspect.getfile(self.__class__)

        ins = sorted(self.inputs(), key=lambda i: i.position if i.position is not None else 0)
        args = ""
        if self.arguments():
            args = " " + " ".join(f"{(a.prefix if a.prefix is not None else '') + ' ' if (a.prefix is not None and a.separate_value_from_prefix) else ''}{a.value}" for a in self.arguments())

        prefixes = " -" + "".join(i.prefix.replace("-", "").replace(" ", "") for i in ins if i.prefix is not None)

        metadata = self.metadata() if self.metadata() else Metadata()
        docker = self.docker()

        base = (self.base_command() if isinstance(self.base_command(), str) else " ".join(self.base_command())) \
            if self.base_command() else ''
        command = base + args + prefixes

        def input_format(t: ToolInput):
            prefix_with_space = ""
            if t.prefix is not None:
                prefix_with_space = (t.prefix + ": ") if (t.separate_value_from_prefix is not False) else t.prefix
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

DOCUMENTATION URL
    {metadata.documentationUrl if metadata.documentationUrl else "No url provided"}

DESCRIPTION
    {metadata.documentation if metadata.documentation else "No documentation provided"}

INPUTS:
    REQUIRED:
{requiredInputs}

    OPTIONAL:
{optionalInputs}

OUTPUTS:
{outputs}
"""
