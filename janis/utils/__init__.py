from typing import Dict, Any, List, Tuple


def first_value(d: Dict):
    return next(iter(d.values()))


def get_value_for_hints_and_ordered_resource_tuple(
    hints: Dict[str, Any], tuples: List[Tuple[str, Dict[str, int]]]
):
    if not hints:
        return None
    for k, d in tuples:
        if k not in hints:
            continue
        v = hints[k]
        if v not in d:
            continue
        return d[v]
    return None


def zip_directory(parent_dir, dir_name):
    import subprocess
    from .logger import Logger
    from os import chdir

    Logger.info("Zipping tools")
    chdir(parent_dir)

    zip_result = subprocess.run(["zip", "-r", f"{dir_name}.zip", f"{dir_name}/"])
    if zip_result.returncode == 0:
        Logger.info("Zipped tools")
    else:
        Logger.critical(zip_result.stderr)
