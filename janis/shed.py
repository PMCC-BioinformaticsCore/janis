from typing import List, Type
import pkg_resources, sys
from inspect import isfunction, ismodule, isabstract, isclass

import janis_core as j

import janis.entrypoints as EP
from janis.registry import TaggedRegistry, Registry


class JanisShed:

    _toolshed = TaggedRegistry("latest")
    _typeshed = Registry()

    _has_been_hydrated = False

    # getters

    @staticmethod
    def get_tool(tool: str, version=None):
        JanisShed.hydrate()
        return JanisShed._toolshed.get(tool, version)

    @staticmethod
    def get_datatype(datatype: str):
        JanisShed.hydrate()
        return JanisShed._typeshed.get(datatype)

    @staticmethod
    def get_all_tools() -> List[List[j.Tool]]:
        JanisShed.hydrate()
        return JanisShed._toolshed.objects()

    @staticmethod
    def get_all_datatypes() -> List[j.DataType]:
        JanisShed.hydrate()
        return JanisShed._typeshed.objects()

    # setters

    @staticmethod
    def add_tool(tool: j.Tool) -> bool:
        v = tool.version()
        if not v:
            t = f"The tool {tool.id()} did not have a version and will not be registered"
            j.Logger.critical(t)
            return False

        return JanisShed._toolshed.register(tool.id(), v, tool)

    @staticmethod
    def add_type(datatype: Type[j.DataType]) -> bool:
        return JanisShed._typeshed.register(datatype.name(), datatype)

    @staticmethod
    def hydrate(force=False, modules: list = None):
        # go get everything
        if JanisShed._has_been_hydrated and not force:
            return

        if not modules:
            modules = []
            # modules.extend(JanisShed._get_datatype_entrypoints())
            modules.extend(JanisShed._get_tool_entrypoints())

        seen_modules = set()
        for m in modules:
            JanisShed.traverse_module(m, seen_modules=seen_modules)

    @staticmethod
    def _get_datatype_entrypoints():
        ep = []
        eps = pkg_resources.iter_entry_points(group=EP.DATATYPES)
        for entrypoint in eps:
            try:
                m = entrypoint.load()
                ep.append(m)
            except ImportError as e:
                t = f"Couldn't import janis data_type extension '{entrypoint.name}': {e}"
                j.Logger.critical(t)
                continue
        return ep

    @staticmethod
    def _get_tool_entrypoints():
        ep = []
        eps = pkg_resources.iter_entry_points(group=EP.TOOLS)
        for entrypoint in eps:
            try:
                m = entrypoint.load()
                ep.append(m)
            except ImportError as e:
                t = f"Couldn't import janis data_type extension '{entrypoint.name}': {e}"
                j.Logger.critical(t)
                continue
        return ep

    @staticmethod
    def traverse_module(module, seen_modules: set):
        if module.__name__ in seen_modules:
            return

        seen_modules.add(module.__name__)

        q = {
            n: cls
            for n, cls in list(module.__dict__.items())
            if not n.startswith("__") and type(cls) != type
        }

        for k in q:
            cls = q[k]
            JanisShed.process_cls(cls, seen_modules)

    @staticmethod
    def process_cls(cls, seen_modules):
        try:
            if ismodule(cls):
                return JanisShed.traverse_module(cls, seen_modules)
            elif isfunction(cls) or isabstract(cls):
                return
            elif not isclass(cls):
                return
            elif issubclass(cls, j.DataType):
                return JanisShed.add_type(cls)
            elif not hasattr(cls, "type") or not callable(cls.type):
                return
            elif cls.type() == j.ToolTypes.Workflow:
                return JanisShed.add_tool(cls())
            elif cls.type() == j.ToolTypes.CommandTool:
                return JanisShed.add_tool(cls())

        except Exception as e:
            j.Logger.log(f"{str(e)} for type {type(cls)}")
