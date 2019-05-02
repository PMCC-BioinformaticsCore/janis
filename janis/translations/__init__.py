from typing import Dict, Optional, Any

from janis.translations.wdl import WdlTranslator
from janis.utils.janisconstants import GITHUB_URL
from janis.translations.exportpath import ExportPathKeywords
from .translationbase import TranslatorBase
from .cwl import CwlTranslator


SupportedTranslation = str


class SupportedTranslations:
    CWL: SupportedTranslation = "cwl"
    WDL: SupportedTranslation = "wdl"


def get_translator(translation: SupportedTranslation) -> TranslatorBase:
    if translation == SupportedTranslations.CWL:
        return CwlTranslator()
    elif translation == SupportedTranslations.WDL:
        return WdlTranslator()

    raise NotImplementedError(f"The requested translation ('{translation}') has not been implemented yet, "
                              f"why not contribute one at '{GITHUB_URL}'.")


def translate_workflow(workflow, translation: SupportedTranslation,
                       to_console=True, with_docker=True, with_resource_overrides=False, to_disk=False,
                       export_path=ExportPathKeywords.default, write_inputs_file=False, should_validate=False,
                       should_zip=True):
    translator = get_translator(translation)
    return translator.translate(
        workflow, to_console=to_console, with_docker=with_docker,
        with_resource_overrides=with_resource_overrides, to_disk=to_disk,
        export_path=export_path, write_inputs_file=write_inputs_file,
        should_validate=should_validate, should_zip=should_zip
    )


def translate_tool(tool, translation: SupportedTranslation, with_docker, with_resource_overrides=False):
    translator = get_translator(translation)
    return translator.stringify_translated_tool(translator.translate_tool(
        tool, with_docker=with_docker,
        with_resource_overrides=with_resource_overrides)
    )


def build_resources_input(workflow, translation: SupportedTranslation, hints) -> str:
    translator = get_translator(translation)
    return translator.stringify_translated_inputs(translator.build_resources_input(workflow, hints))
