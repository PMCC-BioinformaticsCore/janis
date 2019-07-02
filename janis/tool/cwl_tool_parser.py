"""
Generates a CommandTool python stub from a CWL. It's not perfect, but does about 90% of the work for 90% of files.
You will need to import the data_types it uses when copying it into Python, and correct any errors.

It's not parsing requirements, so may not pick up some hints correctly.

PLEASE make sure you adjust your data_types, this only knows about primitives and WON'T HANDLE PRIMITIVES.
"""

CommandToolFormat = """
class {toolClass}(CommandTool):
{required_inputs}

{outputs}

    @staticmethod
    def tool():
        return "{tool}"

    @staticmethod
    def base_command():
        return {baseCommand}

    @staticmethod
    def docker():
        return {docker}

    def doc(self):
        return {doc}

    def arguments(self):
        return [{arguments}]

{optional_inputs}
"""

ToolInputFormat = '{id} = ToolInput("{id}", {typ}{args})'
ToolOutputFormat = '{id} = ToolOutput("{id}", {typ}{args})'
ToolArgumentFormat = 'ToolArgument("{value}"{args})'


def parse_cwl_dict(d):
    tool = d.get("id")

    if tool is None:
        tool = str(input("Tool name: "))

    inputs = []
    optional_inputs = []
    outputs = []
    arguments = []
    docker = d.get("dockerRequirement")
    base_command = d.get("baseCommand")
    doc = d.get("doc")
    if doc is not None:
        doc = doc.replace("\n", "")

    d_inputs = d["inputs"]
    d_outputs = d.get("outputs")
    d_arguments = d.get("arguments")

    if d_inputs is not None:
        for i in d_inputs:
            if isinstance(i, str):
                inp = d_inputs[i]
            else:
                inp = i
                i = inp.get("id")

            args = []
            if isinstance(inp, str):
                typ = inp
            else:

                typ = str(inp.get("type"))

                inpb = inp.get("inputBinding")
                if inpb is not None:
                    if "position" in inpb:
                        args.append("position=" + str(inpb["position"]))
                    if "prefix" in inpb:
                        pre = str(inpb["prefix"])
                        args.append(f'prefix="{pre}"')
                        if pre.endswith("="):
                            args.append("separate_value_from_prefix=False")

                if "doc" in inp:
                    args.append(f'doc="{inp["doc"]}"')

                if "default" in inp:
                    args.append(f'default="{inp["default"]}"')

            optional = "?" in typ or "default" in inp

            if "[]" in typ:
                opt = ", optional=True" if optional else ""
                typ = (
                    "Array("
                    + typ.replace("[]", "").replace("?", "").title()
                    + "()"
                    + opt
                    + ")"
                )
            else:
                opt = "optional=True" if optional else ""
                typ = typ.replace("?", "").title() + "(" + opt + ")"

            s = ToolInputFormat.format(
                id=i, typ=typ, args="".join((", " + x) for x in args)
            )
            if optional:
                optional_inputs.append(s)
            else:
                inputs.append(s)

    if d_outputs is not None:
        for o in d_outputs:
            if isinstance(o, str):
                out = d_outputs[o]
            else:
                out = o
                o = out.get("id")

            typ = out["type"]
            optional = "?" in typ

            if "[]" in typ:
                opt = ", optional=True" if optional else ""
                typ = (
                    "Array("
                    + typ.replace("[]", "").replace("?", "").title()
                    + "()"
                    + opt
                    + ")"
                )
            else:
                opt = "optional=True" if optional else ""
                typ = typ.replace("?", "").title() + "(" + opt + ")"

            args = []

            if "outputBinding" in out:
                if "glob" in out["outputBinding"]:
                    args.append(f"glob='{out['outputBinding']['glob']}'")

            outputs.append(
                ToolOutputFormat.format(
                    id=o, typ=typ, args="".join((", " + x) for x in args)
                )
            )

    if d_arguments is not None:
        for a in d_arguments:
            args = []
            value = a.get("valueFrom")
            if "position" in a:
                args.append("position=" + str(a["position"]))
            if "prefix" in a:
                args.append(f'prefix="{a["prefix"]}"')

            arguments.append(
                ToolArgumentFormat.format(
                    value=value, args="".join((", " + x) for x in args)
                )
            )
    print(
        CommandToolFormat.format(
            toolClass=tool.title(),
            tool=tool,
            baseCommand=base_command,
            required_inputs="\n".join(("    " + i) for i in inputs),
            optional_inputs="\n".join(("    " + i) for i in optional_inputs),
            outputs="\n".join(("    " + o) for o in outputs),
            arguments=", ".join(arguments),
            docker=docker,
            doc=f'"{doc}"',
        )
    )
    # print("\n".join(("    " + i) for i in optional_inputs))


if __name__ == "__main__":
    import argparse, yaml

    parser = argparse.ArgumentParser()
    parser.add_argument("tool")
    args = parser.parse_args()
    print(args)

    with open(args.tool) as tool:
        d = yaml.load(tool)

    parse_cwl_dict(d)
