"""
Welcome to the source code for Janis! This module is mostly a placeholder and
used to bundle the different components of Janis (without crazy circular deps).
It also allows you to add your own tool sheds by including the:

    'janis.extension` entry-point.

Here's an amazing guide on entry-points:
    https://amir.rachum.com/blog/2017/07/28/python-entry-points/

Here are some of the tool sheds that are already built:

    - unix: janis-pipelines.unix (
    - bioinformatics: janis-pipelines.bioinformatics (https://github.com/PMCC-BioinformaticsCore/janis-bioinformatics)

There's also a workflow runner:
    - runner: janis-pipelines.runner (https://github.com/PMCC-BioinformaticsCore/janis-runner)

Janis is a framework creating specialised, simple workflow definitions that are
transpiled to Common Workflow Language or Workflow Definition Language.

_Proudly made on Planet Earth._

"""

import sys, os
import pkg_resources

from janis_core import *

for entrypoint in pkg_resources.iter_entry_points(group="janis.extension"):
    try:
        print(entrypoint.name)
        m = entrypoint.load()
        sys.modules["janis." + entrypoint.name] = m
        globals()[entrypoint.name] = m
    except ImportError as e:
        Logger.critical(f"Couldn't import janis extension '{entrypoint.name}': {e}")
        continue
