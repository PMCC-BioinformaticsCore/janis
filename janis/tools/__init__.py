import sys
from janis_core.utils.logger import Logger
import pkg_resources


for entrypoint in pkg_resources.iter_entry_points(group="janis.tools"):
    try:
        m = entrypoint.load()
        for k, v in m.__dict__.items():
            if k.startswith("_"):
                continue
            globals()[k] = v
    except ImportError as e:
        Logger.critical(f"Couldn't import janis extension '{entrypoint.name}': {e}")
        continue
