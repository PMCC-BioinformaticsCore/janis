from janis_core.utils.logger import Logger
from janis_core.toolbox.entrypoints import DATATYPES
import importlib_metadata

eps = importlib_metadata.entry_points().get(DATATYPES, [])
for entrypoint in eps:
    try:
        m = entrypoint.load()
        for k, v in m.__dict__.items():
            if k.startswith("_"):
                continue
            globals()[k] = v
    except ImportError as e:
        Logger.critical(f"Couldn't import janis types '{entrypoint.name}': {e}")
        continue
