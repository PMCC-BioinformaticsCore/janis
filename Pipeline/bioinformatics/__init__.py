import os
import glob
import importlib

from constants import PROJECT_ROOT_DIR
from Pipeline.types.data_types import DataType
from Pipeline.types.registry import register_type
from Pipeline.tool.commandtool import CommandTool
from Pipeline.tool.registry import register_tool
from Pipeline.utils.logger import Logger

d = os.path.dirname(os.path.abspath(__file__))
for file in glob.glob(os.path.join(d, "**/*.py"), recursive=True):
    # name = os.path.splitext(os.path.basename(file))[0]
    if os.path.basename(file).startswith("__"):
        continue

    name = os.path.splitext(file.replace(PROJECT_ROOT_DIR + "/", ""))[0].replace("/", ".")
    try:
        module = importlib.import_module(name)
        q = {n: cls for n, cls in list(module.__dict__.items()) if not n.startswith("__") and type(cls) != type }
        for cc in q:

            try:
                cls = q[cc]
                if issubclass(cls, CommandTool) and cls != CommandTool:
                    if register_tool(cls):
                        Logger.log("Registered tool: " + cls.tool())

                elif issubclass(cls, DataType) and cls != DataType:
                    if register_type(cls):
                        Logger.log("Registered type: " + cls.name())
            except:
                continue

    except:
        continue
