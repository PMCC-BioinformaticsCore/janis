from janis.utils.janisconstants import GITHUB_URL
from . import cwl, wdl
from .wdl import dump_wdl as wdl_dump_workflow

SupportedTranslation = str


class SupportedTranslations:
    CWL: SupportedTranslation = "cwl"
    WDL: SupportedTranslation = "wdl"


def dump_translation(workflow, translation: SupportedTranslation, to_console=True, to_disk=False, with_docker=True,
                     with_hints=False, with_resource_overrides=False):
    if translation == SupportedTranslations.CWL:
        return cwl.dump_cwl(workflow, to_console=to_console, to_disk=to_disk, with_docker=with_docker,
                                 with_hints=with_hints, with_resource_overrides=with_resource_overrides)

    elif translation == SupportedTranslations.WDL:
        return wdl.dump_wdl(workflow, to_console=to_console, to_disk=to_disk, with_docker=with_docker,
                                 with_hints=with_hints, with_resource_overrides=with_resource_overrides)

    else:
        raise NotImplementedError(f"The requested translation ('{translation}') has not been implemented yet, "
                                  f"why not contribute one at '{GITHUB_URL}'.")


def translate_tool(tool, translation: SupportedTranslation, with_docker):
    if translation == SupportedTranslations.CWL:
        return cwl.translate_tool_str(tool, with_docker=with_docker)

    elif translation == SupportedTranslations.WDL:
        return wdl.translate_tool(tool, with_docker=with_docker)

    else:
        raise NotImplementedError(f"The requested translation ('{translation}') has not been implemented yet, "
                                  f"why not contribute one at '{GITHUB_URL}'.")
