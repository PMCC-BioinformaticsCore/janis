import os
import glob
import importlib

from types.data_types import DataType
from Tool.tool import Tool
from pipeline_definition.types.type_registry import register_tool, register_type
from constants import PROJECT_ROOT_DIR
from utils.logger import Logger

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
                if issubclass(cls, Tool) and cls != Tool:
                    if register_tool(cls):
                        Logger.log("Registered tool: " + cls.tool())

                elif issubclass(cls, DataType) and cls != DataType:
                    if register_type(cls):
                        Logger.log("Registered type: " + cls.name())
            except:
                continue

    except:
        continue
