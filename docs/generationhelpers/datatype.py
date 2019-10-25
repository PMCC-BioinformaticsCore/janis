from datetime import date
from janis_core import DataType


def prepare_data_type(dt: DataType):
    dt_name = dt.name()
    secondary = ""

    if dt.secondary_files():
        secondary = "Secondary files: " + ", ".join(
            f"``{s}``" for s in dt.secondary_files()
        )

    return f"""
{dt_name}
{"=" * len(dt_name)}

{secondary}

Documentation
-------------

{dt.doc()}

*This page was automatically generated on {date.today().strftime("%Y-%m-%d")}*.
"""
