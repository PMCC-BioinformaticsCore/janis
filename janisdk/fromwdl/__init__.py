def add_fromwdl_args(parser):
    parser.description = (
        "Parse a WDL command line tool / workflow and write it as a Janis file"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Directory to output the workflow / tools to, otherwise this is written to stdout",
    )

    parser.add_argument("wdlfile", help="The path to the WDL file")

    parser.add_argument(
        "translation", default="janis", choices=["cwl", "wdl", "janis"], nargs="?"
    )

    return parser


def do_fromwdl(args):
    from janis_core import WdlParser, Logger


    Logger.info(f"Loading WDL file: {args.wdlfile}")
    tool = WdlParser.from_doc(args.wdlfile)

    Logger.info(f"Loaded {tool.type()}: {tool.versioned_id()}")

    translated = tool.translate(
        args.translation,
        to_console=args.output is None,
        to_disk=args.output is not None,
        export_path=args.output,
    )

    return translated
