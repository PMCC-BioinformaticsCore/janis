import os
import glob
import importlib
from inspect import isclass, isfunction, signature, isabstract

from Pipeline import Workflow
from Pipeline.tool.tool import Tool
from constants import PROJECT_ROOT_DIR
from Pipeline.types.data_types import DataType
from Pipeline.types.registry import register_type
from Pipeline.tool.commandtool import CommandTool
from Pipeline.tool.registry import register_tool, get_tools
from Pipeline.utils.logger import Logger

def try_register_type(cls):
    try:
        if isabstract(cls): return
        # if isfunction(cls) and len(signature(cls).parameters) == 0:
        #     Logger.log(str(cls) + " is an empty-param function")
        #     Logger.mute()
        #     resolved = cls()
        #     Logger.unmute()
        #     if isinstance(resolved, Workflow):
        #         if register_tool(cls, resolved.id()):
        #             Logger.log("Registered wf as tool: " + resolved.id())

        elif isclass(cls) and issubclass(cls, CommandTool) and cls != CommandTool:
            Logger.log("attempting to register " + str(cls))
            if register_tool(cls):
                Logger.log("Registered tool: " + cls.tool())

        elif isclass(cls) and issubclass(cls, Workflow) and cls != Workflow:
            identifier = cls().id()
            if register_tool(cls, identifier):
                Logger.log("Registered wf as tool: " + identifier)

        elif isclass(cls) and issubclass(cls, DataType) and cls != DataType:
            if register_type(cls):
                Logger.log("Registered type: " + cls.name())

    except Exception as e:
        import traceback
        Logger.log_ex(e)
        Logger.warn(traceback.format_exc())

ignore_files = set(["regeneratedocumentation.py"])

d = os.path.dirname(os.path.abspath(__file__))
Logger.log("Locating modules from " + d)
files = list(glob.glob(os.path.join(d, "**/*.py"), recursive=True))
Logger.info(f"Finding modules in {len(files)} files")
for file in files:
    if os.path.basename(file).startswith("__"):
        continue
    if os.path.basename(file) in ignore_files:
        continue


    name = os.path.splitext(file.replace(PROJECT_ROOT_DIR + "/", ""))[0].replace("/", ".")
    try:
        module = importlib.import_module(name)
        q = {n: cls for n, cls in list(module.__dict__.items()) if not n.startswith("__") and type(cls) != type}
        for cc in q:
            try_register_type(q[cc])

    except:
        continue
