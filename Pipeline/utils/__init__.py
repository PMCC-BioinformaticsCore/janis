from typing import Dict


def first_value(d: Dict):
    return next(iter(d.values()))