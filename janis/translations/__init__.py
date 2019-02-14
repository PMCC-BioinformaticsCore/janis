from janis.utils.janisconstants import GITHUB_URL
from .cwl import dump_cwl as cwl_dump_workflow
from .wdl import dump_wdl as wdl_dump_workflow

SupportedTranslation = str

class SupportedTranslations:
    CWL: SupportedTranslation = "cwl"
    WDL: SupportedTranslation = "wdl"


def dump_translation(workflow, translation: SupportedTranslation, to_console=True, to_disk=False, with_docker=True,
                     with_hints=False, with_resource_overrides=False):
    if translation == SupportedTranslations.CWL:
        return cwl_dump_workflow(workflow, to_console=to_console, to_disk=to_disk, with_docker=with_docker,
                                 with_hints=with_hints, with_resource_overrides=with_resource_overrides)

    elif translation == SupportedTranslations.WDL:
        return wdl_dump_workflow(workflow, to_console=to_console, to_disk=to_disk, with_docker=with_docker,
                                 with_hints=with_hints, with_resource_overrides=with_resource_overrides)

    else:
        raise NotImplementedError(f"The requested translation ('{translation}') has not been implemented yet, "
                                  f"why not contribute one at '{GITHUB_URL}'.")
