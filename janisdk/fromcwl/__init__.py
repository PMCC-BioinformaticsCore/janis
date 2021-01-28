def add_fromcwl_args(parser):
    parser.description = (
        "Parse a CWL command line tool / workflow and write it as a Janis file"
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Directory to output the workflow / tools to, otherwise this is written to stdout",
    )

    parser.add_argument("cwlfile", help="The path to the CWL file")

    parser.add_argument(
        "translation", default="janis", choices=["cwl", "wdl", "janis"], nargs="?"
    )

    return parser


def do_fromcwl(args):
    from janis_core import CWlParser, Logger

    Logger.info(f"Loading CWL file: {args.cwlfile}")
    tool = CWlParser.from_doc(args.cwlfile)

    Logger.info(f"Loaded {tool.type()}: {tool.versioned_id()}")

    translated = tool.translate(
        args.translation,
        to_console=args.output is None,
        to_disk=args.output is not None,
        export_path=args.output,
    )

    return translated
