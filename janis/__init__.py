"""
Welcome to the source code for Janis!

Janis is a framework creating specialised, simple workflow definitions that are
transpiled to Common Workflow Language or Workflow Definition Language.

Below you'll find the classes you can access by importing Janis:

>>  import janis as j

Some noteworthy classes are Workflow, CommandTool and the janis.translations module

Some terminology:
    - Edge:
        - An edge may exist between 2 nodes, it represents a dependency
        - Every edge will have a source_map which is indexed by the tag of the input on the finish node
        - This source_map value is singular or a list, as you can connect multiple sources to one input



_Proudly made on Planet Earth._

"""

import sys, os
import pkg_resources

from janis_core import *

runner = None

for entrypoint in pkg_resources.iter_entry_points(group="janis.extension"):
    try:
        print(entrypoint.name)
        m = entrypoint.load()
        sys.modules["janis." + entrypoint.name] = m
        globals()[entrypoint.name] = m
    except ImportError as e:
        Logger.critical(f"Couldn't import janis extension '{entrypoint.name}': {e}")
        continue
