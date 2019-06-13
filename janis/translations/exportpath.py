class ExportPathKeywords:
    workflow_spec = "{language}"
    workflow_name = "{name}"

    default = f"~/Desktop/{workflow_name}/{workflow_spec}/"
    default_no_spec = f"~/Desktop/{workflow_name}/"

    @staticmethod
    def resolve(path, workflow_spec, workflow_name):
        from os.path import expanduser

        if ExportPathKeywords.workflow_spec in path and workflow_spec is None:
            raise Exception(f"path ('{path}') contained parameter {ExportPathKeywords.workflow_spec} "
                            "but caller of .resolve did not pass language")

        if ExportPathKeywords.workflow_name in path and workflow_name is None:
            raise Exception(f"path ('{path}') contained parameter {ExportPathKeywords.workflow_name} "
                            "but caller of .resolve did not pass workflow name")

        return expanduser(path) \
            .replace(ExportPathKeywords.workflow_spec, workflow_spec.lower() if workflow_spec else "") \
            .replace(ExportPathKeywords.workflow_name, workflow_name.lower() if workflow_name else "")
