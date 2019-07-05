from enum import Enum
from typing import Dict

from janis.translations.exportpath import ExportPathKeywords
from janis.translations.wdl import WdlTranslator
from janis.__meta__ import GITHUB_URL
from .cwl import CwlTranslator
from .translationbase import TranslatorBase

SupportedTranslation = str


class SupportedTranslations(Enum):
    CWL = "cwl"
    WDL = "wdl"

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return str(self) == str(other)


def get_translator(translation: SupportedTranslation) -> TranslatorBase:
    if translation == SupportedTranslations.CWL:
        return CwlTranslator()
    elif translation == SupportedTranslations.WDL:
        return WdlTranslator()

    raise NotImplementedError(
        f"The requested translation ('{translation}') has not been implemented yet, "
        f"why not contribute one at '{GITHUB_URL}'."
    )


def translate_workflow(
    workflow,
    translation: SupportedTranslation,
    to_console=True,
    with_docker=True,
    with_resource_overrides=False,
    to_disk=False,
    write_inputs_file=True,
    export_path=ExportPathKeywords.default,
    should_validate=False,
    should_zip=True,
    merge_resources=False,
    hints=None,
    allow_null_if_not_optional=True,
    additional_inputs: Dict = None,
):
    translator = get_translator(translation)
    return translator.translate(
        workflow,
        to_console=to_console,
        with_docker=with_docker,
        with_resource_overrides=with_resource_overrides,
        to_disk=to_disk,
        export_path=export_path,
        write_inputs_file=write_inputs_file,
        should_validate=should_validate,
        should_zip=should_zip,
        merge_resources=merge_resources,
        hints=hints,
        allow_null_if_not_optional=allow_null_if_not_optional,
        additional_inputs=additional_inputs,
    )


def translate_tool(
    tool,
    translation: SupportedTranslation,
    to_console=True,
    with_docker=True,
    with_resource_overrides=False,
) -> str:
    translator = get_translator(translation)
    tool_out = translator.stringify_translated_tool(
        translator.translate_tool(
            tool,
            with_docker=with_docker,
            with_resource_overrides=with_resource_overrides,
        )
    )

    if to_console:
        print(tool_out)

    return tool_out


def build_resources_input(workflow, translation: SupportedTranslation, hints) -> str:
    translator = get_translator(translation)
    return translator.stringify_translated_inputs(
        translator.build_resources_input(workflow, hints)
    )
