from janis.utils.logger import Logger
from os import path as p, getcwd


class ExportPathKeywords:
    workflow_spec = "{language}"
    workflow_name = "{name}"

    default = "./{name}"  # f"~/Desktop/{workflow_name}/{workflow_spec}/"
    default_no_spec = "./{name}"  # f"~/Desktop/{workflow_name}/"

    @staticmethod
    def resolve(path, workflow_spec, workflow_name):

        if not path:
            Logger.critical("Output path was invalid, changed to working directory")
            path = "."

        if ExportPathKeywords.workflow_spec in path and workflow_spec is None:
            raise Exception(
                f"path ('{path}') contained parameter {ExportPathKeywords.workflow_spec} "
                "but caller of .resolve did not pass language"
            )

        if ExportPathKeywords.workflow_name in path and workflow_name is None:
            raise Exception(
                f"path ('{path}') contained parameter {ExportPathKeywords.workflow_name} "
                "but caller of .resolve did not pass workflow name"
            )

        if (len(path) == 1 and path == ".") or path[:2] == "./":
            path = getcwd() + (path[1:] if len(path) > 1 else "")

        return (
            p.expanduser(path)
            .replace(
                ExportPathKeywords.workflow_spec,
                workflow_spec.lower() if workflow_spec else "",
            )
            .replace(
                ExportPathKeywords.workflow_name,
                workflow_name.lower() if workflow_name else "",
            )
        )
