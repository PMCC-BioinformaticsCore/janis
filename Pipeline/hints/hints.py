# process for adding new hints!
from typing import Union

recognisedHints = {
    "captureType": ["targeted", "exome", "30x", "90x", "300x"]
}


def get_cwl_schema_for_recognised_hints():
    import cwlgen.cwlgen as cwl

    schema = cwl.CommandInputRecordSchema("hints")

    def prepare_hint(hint, value):
        name = hint + "_schema"

        if isinstance(value, list):
            # assume is an enum
            return cwl.CommandInputEnumSchema(label=hint, name=name, symbols=value)
        elif value == "array":
            return cwl.CommandInputArraySchema("string", label=hint)
        else:
            return "string"

    schema.fields = [cwl.CommandInputRecordSchema.CommandInputRecordField(hint, ["null", prepare_hint(hint, recognisedHints[hint])])
        for hint in recognisedHints]

    return schema
